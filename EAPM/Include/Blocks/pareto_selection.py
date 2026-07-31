"""
Module containing the Pareto selection block for the EAPM plugin.

Reproduces the UPO notebook's multi-objective ranking and final selection
(cells 55-67): the per-model metrics are merged into one table, every model gets
a Pareto rank, and exactly N models are picked by rank with crowding distance
breaking ties inside the last front.

The notebook's objectives are ``Q_Average`` (maximised) and the per-residue
Rosetta score (minimised), but any numeric columns can be used, in any number.

Both the ranking (cell 56) and the crowding-distance selection (cell 63) are
pure numpy, so this block has no external tool dependency.
"""

from HorusAPI import (
    Extensions,
    PluginBlock,
    PluginVariable,
    VariableList,
    VariableTypes,
)

# ==========================#
# Variable inputs
# ==========================#
metricsFile = PluginVariable(
    id="metrics_file",
    name="Metrics file",
    description="CSV of per-model metrics indexed by model (e.g. 'bioemu_metrics.csv').",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
extraMetricsFile = PluginVariable(
    id="extra_metrics_file",
    name="Extra metrics file (optional)",
    description="A second per-model CSV merged on the model name "
    "(e.g. 'rosetta_score_summary.csv', which carries score_per_residue and length).",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["csv"],
)
originFile = PluginVariable(
    id="origin_file",
    name="Model origin (optional)",
    description="JSON {model: origin} from 'Combine Selected Sequences'. Added as a "
    "Model_Origin column for reference.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
objectiveColumnVariable = PluginVariable(
    id="objective_column",
    name="Column",
    description="Name of the metric column to optimise.",
    type=VariableTypes.STRING,
)
objectiveMaximizeVariable = PluginVariable(
    id="objective_maximize",
    name="Maximize",
    description="Enable to maximise this objective; leave off to minimise it.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
objectivesList = VariableList(
    id="objectives",
    name="Objectives",
    description="Metrics to optimise. Defaults to the notebook's pair: Q_Average "
    "maximised and score_per_residue minimised.",
    prototypes=[objectiveColumnVariable, objectiveMaximizeVariable],
    defaultValue=[
        {"objective_column": "Q_Average", "objective_maximize": True},
        {"objective_column": "score_per_residue", "objective_maximize": False},
    ],
)
numSelectVariable = PluginVariable(
    id="n_select",
    name="Models to select",
    description="How many models to select by Pareto rank and crowding distance.",
    type=VariableTypes.INTEGER,
    defaultValue=10,
)

# ==========================#
# Variable outputs
# ==========================#
rankedFile = PluginVariable(
    id="ranked_file",
    name="Ranked models",
    description="CSV with every model, its Pareto_Rank and whether it was selected.",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
selectedFile = PluginVariable(
    id="selected_file",
    name="Selected models table",
    description="CSV with the selected models and all their metrics.",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
selectedModelsFile = PluginVariable(
    id="selected_models_file",
    name="Selected models",
    description="JSON list with the names of the selected models.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
plotFile = PluginVariable(
    id="plot_file",
    name="Pareto plot",
    description="PNG comparing the objectives against the other metrics, coloured by "
    "Pareto rank with the selected models highlighted.",
    type=VariableTypes.FILE,
    allowedValues=["png"],
)


def _pareto_ranks(objectives):
    """
    Assign a Pareto rank to every row of ``objectives`` (all columns minimised).

    Rank 1 is the non-dominated front, rank 2 the front once rank 1 is removed,
    and so on. Mirrors the notebook's compute_pareto_ranks.
    """
    import numpy as np

    count = objectives.shape[0]
    ranks = np.zeros(count, dtype=int)
    remaining = set(range(count))
    rank = 1

    while remaining:
        front = []
        for i in remaining:
            dominated = False
            for j in remaining:
                if i == j:
                    continue
                if np.all(objectives[j] <= objectives[i]) and np.any(
                    objectives[j] < objectives[i]
                ):
                    dominated = True
                    break
            if not dominated:
                front.append(i)

        if not front:  # defensive: never loop forever on degenerate input
            front = list(remaining)

        for i in front:
            ranks[i] = rank
        remaining -= set(front)
        rank += 1

    return ranks


def _crowding_distance(objectives):
    """Crowding distance of each point in objective space (all minimised)."""
    import numpy as np

    count, n_obj = objectives.shape
    distances = np.zeros(count)
    if count <= 2:
        return np.full(count, np.inf)

    for j in range(n_obj):
        order = np.argsort(objectives[:, j])
        distances[order[0]] = distances[order[-1]] = np.inf
        column = objectives[order, j]
        spread = column[-1] - column[0]
        if spread == 0:
            continue
        for i in range(1, count - 1):
            distances[order[i]] += (column[i + 1] - column[i - 1]) / spread

    return distances


def _select_by_crowding(ranks, objectives, n_select):
    """Pick ``n_select`` indices by Pareto rank, using crowding distance within a front."""
    import numpy as np

    selected = []
    for rank in sorted(set(ranks.tolist())):
        front = np.where(ranks == rank)[0].tolist()
        slots = n_select - len(selected)
        if slots <= 0:
            break
        if len(front) <= slots:
            selected.extend(front)
        else:
            distances = _crowding_distance(objectives[front])
            best = np.argsort(-distances)[:slots]
            selected.extend([front[i] for i in best])
    return selected


def _build_plot(frame, objective_columns, output_path):
    """Scatter the first objective against every other metric, coloured by rank."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import numpy as np

    y_column = objective_columns[0]
    numeric = [
        c
        for c in frame.columns
        if c != y_column
        and c not in ("Pareto_Rank", "Selected")
        and np.issubdtype(frame[c].dtype, np.number)
    ]
    if not numeric:
        numeric = [y_column]

    ranks = frame["Pareto_Rank"].to_numpy()
    norm = mcolors.Normalize(vmin=ranks.min(), vmax=ranks.max())
    cmap = plt.get_cmap("viridis")

    cols = min(3, len(numeric))
    rows = (len(numeric) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)

    selected_mask = frame["Selected"].to_numpy(dtype=bool)

    for index, column in enumerate(numeric):
        axis = axes[index // cols][index % cols]
        axis.scatter(
            frame.loc[~selected_mask, column],
            frame.loc[~selected_mask, y_column],
            c=[cmap(norm(r)) for r in ranks[~selected_mask]],
            s=45,
        )
        axis.scatter(
            frame.loc[selected_mask, column],
            frame.loc[selected_mask, y_column],
            facecolors="none",
            edgecolors="red",
            linewidths=1.8,
            s=130,
            label="selected",
        )
        axis.set_xlabel(column)
        axis.set_ylabel(y_column)
        axis.set_title(f"{y_column} vs {column}", fontsize=10)
        axis.grid(True, alpha=0.3)

    for index in range(len(numeric), rows * cols):
        fig.delaxes(axes[index // cols][index % cols])

    handles, labels = axes[0][0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right")

    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    fig.colorbar(mappable, ax=axes.ravel().tolist(), label="Pareto rank", shrink=0.6)

    fig.suptitle("Pareto ranking (red = selected)", fontsize=14)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def pareto_selection(block: PluginBlock):
    """
    Rank the models by Pareto dominance and select the best N.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import base64
    import json
    import os

    import numpy as np
    import pandas as pd

    # pylint: enable=import-outside-toplevel

    metrics_path = block.inputs.get(metricsFile.id, None)
    if not metrics_path or not os.path.isfile(metrics_path):
        raise ValueError("A valid metrics CSV must be provided.")

    frame = pd.read_csv(metrics_path, index_col=0)
    frame.index = frame.index.astype(str)

    extra_path = block.inputs.get(extraMetricsFile.id, None)
    if extra_path and os.path.isfile(extra_path):
        extra = pd.read_csv(extra_path, index_col=0)
        extra.index = extra.index.astype(str)
        overlapping = [c for c in extra.columns if c in frame.columns]
        if overlapping:
            print(f"Dropping duplicated columns from the extra metrics: {overlapping}")
            extra = extra.drop(columns=overlapping)
        frame = frame.join(extra, how="outer")
        print(f"Merged extra metrics: {list(extra.columns)}")

    origin_path = block.inputs.get(originFile.id, None)
    if origin_path and os.path.isfile(origin_path):
        with open(origin_path) as jf:
            origins = json.load(jf)
        frame["Model_Origin"] = [origins.get(str(m), "unknown") for m in frame.index]

    # ---- Objectives ----
    raw_objectives = block.variables.get(objectivesList.id) or []
    objectives = []
    for entry in raw_objectives:
        column = (entry.get("objective_column") or "").strip()
        if not column:
            continue
        objectives.append((column, bool(entry.get("objective_maximize", False))))

    if not objectives:
        raise ValueError("At least one objective must be defined.")

    missing = [c for c, _ in objectives if c not in frame.columns]
    if missing:
        raise ValueError(
            f"Objective column(s) {missing} are not in the metrics table. "
            f"Available columns: {list(frame.columns)}"
        )

    objective_columns = [c for c, _ in objectives]
    maximize = [m for _, m in objectives]

    # Rows with a missing objective cannot be compared; drop them explicitly.
    usable = frame.dropna(subset=objective_columns)
    dropped = len(frame) - len(usable)
    if dropped:
        print(f"Warning: {dropped} model(s) dropped for missing objective values.")
    if usable.empty:
        raise ValueError("No models have values for every objective.")

    n_select = int(block.variables.get(numSelectVariable.id, 10) or 10)
    if n_select > len(usable):
        print(f"Only {len(usable)} model(s) available; selecting all of them.")
        n_select = len(usable)

    # Negate maximised objectives so everything is minimised
    data = usable[objective_columns].to_numpy(dtype=float)
    for index, do_max in enumerate(maximize):
        if do_max:
            data[:, index] = -data[:, index]

    ranks = _pareto_ranks(data)
    selected_positions = _select_by_crowding(ranks, data, n_select)

    result = usable.copy()
    result["Pareto_Rank"] = ranks
    result["Selected"] = False
    result.iloc[
        selected_positions, result.columns.get_loc("Selected")
    ] = True

    result = result.sort_values(
        ["Selected", "Pareto_Rank"], ascending=[False, True]
    )

    selected = result[result["Selected"]]

    described = ", ".join(
        f"{c} ({'max' if m else 'min'})" for c, m in objectives
    )
    print(f"Objectives: {described}")
    print(
        f"Ranked {len(result)} model(s) into {int(ranks.max())} front(s); "
        f"selected {len(selected)}."
    )

    ranked_output = "pareto_ranked_models.csv"
    result.to_csv(ranked_output)

    selected_output = "pareto_selected_models.csv"
    selected.to_csv(selected_output)

    selected_models_output = "pareto_selected_models.json"
    with open(selected_models_output, "w") as jf:
        json.dump([str(m) for m in selected.index], jf, indent=2)

    plot_output = "pareto_selection_plot.png"
    _build_plot(result, objective_columns, plot_output)

    with open(plot_output, "rb") as img:
        encoded = base64.b64encode(img.read()).decode("ascii")

    html = f"""
    <style>
      .eapm-wrap {{ font-family: sans-serif; }}
      .eapm-scroll {{ max-height: 420px; overflow: auto; margin-bottom: 24px; }}
      table.eapm-table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
      table.eapm-table th, table.eapm-table td {{
          border: 1px solid #ddd; padding: 4px 8px; text-align: right;
      }}
      table.eapm-table th {{ background: #f2f2f2; position: sticky; top: 0; }}
    </style>
    <div class="eapm-wrap">
      <h3>Pareto selection</h3>
      <p>Objectives: <b>{described}</b> &middot; {len(result)} models ranked into
         {int(ranks.max())} fronts &middot; <b>{len(selected)}</b> selected</p>
      <div style="overflow-x:auto; margin-bottom:24px">
        <img src="data:image/png;base64,{encoded}" style="max-width:100%"/>
      </div>
      <h3>Selected models</h3>
      <div class="eapm-scroll">
        {selected.round(3).to_html(classes="eapm-table", border=0)}
      </div>
      <h3>All models</h3>
      <div class="eapm-scroll">
        {result.round(3).to_html(classes="eapm-table", border=0)}
      </div>
    </div>
    """
    Extensions().loadHTML(html, title="Pareto selection")

    block.setOutput(rankedFile.id, ranked_output)
    block.setOutput(selectedFile.id, selected_output)
    block.setOutput(selectedModelsFile.id, selected_models_output)
    block.setOutput(plotFile.id, plot_output)


paretoSelectionBlock = PluginBlock(
    category="Model Ranking",
    name="Pareto Selection",
    id="pareto_selection",
    description="Rank models by Pareto dominance over any number of objectives and "
    "select the best N using crowding distance.",
    inputs=[metricsFile, extraMetricsFile, originFile],
    variables=[objectivesList, numSelectVariable],
    outputs=[rankedFile, selectedFile, selectedModelsFile, plotFile],
    action=pareto_selection,
)
