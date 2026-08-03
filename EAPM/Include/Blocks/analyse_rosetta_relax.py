"""
Module containing the Rosetta relax analysis block for the EAPM plugin.

Reproduces the UPO notebook's relax analysis (cells 30-34)::

    jobs = selected_models.analyseRosettaCalculation('relax_selected_models',
                                                      return_jobs=True,
                                                      skip_finished=True)
    if jobs:
        bsc_calculations.local.parallel(jobs, script_name='relax_commands')
    relax_data = selected_models.analyseRosettaCalculation('relax_selected_models')

The extraction runs in two phases: ``return_jobs=True`` yields the per-model
commands that mine the Rosetta silent files, those are run locally in parallel,
and a second call collects the results into a DataFrame.

Note that ``bsc_calculations.local.parallel`` only *writes* the scripts, and the
main script it generates launches the numbered ones with ``nohup ... &`` without
waiting. This block therefore runs the numbered scripts itself and waits for all
of them, so the DataFrame is only read once every job has finished.
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
modelsFolder = PluginVariable(
    id="models_folder",
    name="Models folder",
    description="Folder with the models that were relaxed (e.g. 'selected_models').",
    type=VariableTypes.FOLDER,
)
relaxFolder = PluginVariable(
    id="relax_folder",
    name="Relax folder",
    description="Rosetta relax output folder produced by the 'Rosetta Relax' block.",
    type=VariableTypes.FOLDER,
)
lengthsFile = PluginVariable(
    id="lengths_file",
    name="Sequence lengths (optional)",
    description="Optional JSON {model: length} (e.g. 'sequence_lengths.json' from the "
    "'Sequence Length Distribution' block). When given, a score-per-residue panel "
    "is added to the plot and the normalised score is written to the summary.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
cpusVariable = PluginVariable(
    id="cpus",
    name="CPUs",
    description="How many extraction jobs to run in parallel locally.",
    type=VariableTypes.INTEGER,
    defaultValue=4,
)
skipFinishedVariable = PluginVariable(
    id="skip_finished",
    name="Skip finished",
    description="Do not regenerate extraction jobs whose output already exists.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)
overwriteVariable = PluginVariable(
    id="overwrite",
    name="Overwrite",
    description="Force re-extraction of the Rosetta data.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
energyByResidueVariable = PluginVariable(
    id="energy_by_residue",
    name="Energy by residue",
    description="Also extract per-residue energies (slower).",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
verboseVariable = PluginVariable(
    id="verbose",
    name="Verbose",
    description="Show the extraction output in the logs.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
maxRowsVariable = PluginVariable(
    id="max_display_rows",
    name="Rows per page",
    description="Number of rows shown per page in the 'Relax data' table. All the "
    "rows are available in the view, this only sets the page size.",
    type=VariableTypes.INTEGER,
    defaultValue=200,
)
plotTypeVariable = PluginVariable(
    id="plot_type",
    name="Plot type",
    description="How to show each model's score distribution.",
    type=VariableTypes.STRING,
    defaultValue="violin",
    allowedValues=["violin", "box", "strip"],
)
maxModelsPlotVariable = PluginVariable(
    id="max_models_plot",
    name="Max models in plot",
    description="Only the best-scoring N models are drawn, keeping the x-axis "
    "readable. Set to 0 to plot every model.",
    type=VariableTypes.INTEGER,
    defaultValue=60,
)

# ==========================#
# Variable outputs
# ==========================#
relaxDataFile = PluginVariable(
    id="relax_data_file",
    name="Relax data",
    description="CSV with the full Rosetta relax data (one row per model/pose).",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
scoreSummaryFile = PluginVariable(
    id="score_summary_file",
    name="Score summary",
    description="CSV with per-model score statistics (mean, min, max, std, count).",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
meanScoresFile = PluginVariable(
    id="mean_scores_file",
    name="Mean scores",
    description="JSON mapping each model to its average Rosetta score.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
plotFile = PluginVariable(
    id="plot_file",
    name="Score plot",
    description="PNG with the per-model Rosetta score distributions.",
    type=VariableTypes.FILE,
    allowedValues=["png"],
)


_BASE_CSS = """
  .eapm-wrap { font-family: sans-serif; margin: 0 auto; padding: 12px 16px; }
  .eapm-wrap h3 { margin-bottom: 4px; }
  .eapm-scroll { max-height: 420px; overflow: auto; margin-bottom: 24px; }
  table.eapm-table { border-collapse: collapse; width: 100%; font-size: 12px; }
  table.eapm-table th, table.eapm-table td {
      border: 1px solid #ddd; padding: 4px 8px; text-align: right;
  }
  table.eapm-table th { background: #f2f2f2; position: sticky; top: 0; }
"""

_TABLE_CSS = """
  .eapm-controls {
      display: flex; gap: 12px; align-items: center; flex-wrap: wrap;
      margin-bottom: 8px;
  }
  .eapm-controls input, .eapm-controls select {
      font-size: 13px; padding: 4px 6px; border: 1px solid #ccc; border-radius: 4px;
  }
  .eapm-controls input { min-width: 220px; }
  .eapm-count { color: #666; font-size: 12px; }
  .eapm-pager { display: flex; gap: 6px; align-items: center; margin-top: 8px; }
  .eapm-pager button {
      font-size: 13px; padding: 3px 10px; border: 1px solid #ccc;
      border-radius: 4px; background: #fafafa; cursor: pointer;
  }
  .eapm-pager button:disabled { opacity: 0.45; cursor: default; }
  table.eapm-table th.eapm-sortable { cursor: pointer; user-select: none; }
  table.eapm-table th.eapm-sortable:hover { background: #e6e6e6; }
  table.eapm-table tbody tr:nth-child(even) { background: #fafafa; }
"""

# Rendered client side so the whole dataset stays usable: filter, sort and page
# through every row without regenerating the page.
_TABLE_JS = """
  const rows = PAYLOAD.data;
  const columns = PAYLOAD.columns;

  let filtered = rows;
  let sortCol = -1;
  let sortAsc = true;
  let page = 0;
  let pageSize = PAGE_SIZE;

  const search = document.getElementById('eapm-search');
  const pageSizeSelect = document.getElementById('eapm-page-size');
  const head = document.getElementById('eapm-head');
  const body = document.getElementById('eapm-body');
  const count = document.getElementById('eapm-count');
  const pageInfo = document.getElementById('eapm-page-info');
  const prev = document.getElementById('eapm-prev');
  const next = document.getElementById('eapm-next');

  function renderHead() {
    const tr = document.createElement('tr');
    columns.forEach((name, index) => {
      const th = document.createElement('th');
      th.className = 'eapm-sortable';
      th.textContent = name;
      th.addEventListener('click', () => {
        if (sortCol === index) {
          sortAsc = !sortAsc;
        } else {
          sortCol = index;
          sortAsc = true;
        }
        applySort();
        page = 0;
        render();
      });
      tr.appendChild(th);
    });
    head.appendChild(tr);
  }

  function updateHeadLabels() {
    Array.from(head.querySelectorAll('th')).forEach((th, index) => {
      const arrow = index === sortCol ? (sortAsc ? ' \\u25b2' : ' \\u25bc') : '';
      th.textContent = columns[index] + arrow;
    });
  }

  function applySort() {
    if (sortCol < 0) {
      return;
    }
    const direction = sortAsc ? 1 : -1;
    filtered = filtered.slice().sort((a, b) => {
      const left = a[sortCol];
      const right = b[sortCol];
      if (left === null || left === undefined) return 1;
      if (right === null || right === undefined) return -1;
      if (typeof left === 'number' && typeof right === 'number') {
        return (left - right) * direction;
      }
      return String(left).localeCompare(String(right), undefined, {
        numeric: true
      }) * direction;
    });
  }

  function applyFilter() {
    const term = search.value.trim().toLowerCase();
    filtered = !term
      ? rows
      : rows.filter((row) =>
          row.some((cell) => cell !== null && String(cell).toLowerCase().includes(term))
        );
    applySort();
  }

  function render() {
    const total = filtered.length;
    const pages = Math.max(1, Math.ceil(total / pageSize));
    page = Math.min(page, pages - 1);

    const start = page * pageSize;
    const slice = filtered.slice(start, start + pageSize);

    body.textContent = '';
    const fragment = document.createDocumentFragment();
    slice.forEach((row) => {
      const tr = document.createElement('tr');
      row.forEach((cell) => {
        const td = document.createElement('td');
        td.textContent = cell === null || cell === undefined ? '' : cell;
        tr.appendChild(td);
      });
      fragment.appendChild(tr);
    });
    body.appendChild(fragment);

    count.textContent = total === rows.length
      ? total + ' rows'
      : total + ' of ' + rows.length + ' rows';
    pageInfo.textContent = total
      ? 'Rows ' + (start + 1) + '-' + Math.min(start + pageSize, total) +
        ' (page ' + (page + 1) + ' of ' + pages + ')'
      : 'No matching rows';
    prev.disabled = page === 0;
    next.disabled = page >= pages - 1;
    updateHeadLabels();
  }

  search.addEventListener('input', () => {
    applyFilter();
    page = 0;
    render();
  });

  pageSizeSelect.addEventListener('change', () => {
    pageSize = parseInt(pageSizeSelect.value, 10);
    page = 0;
    render();
  });

  prev.addEventListener('click', () => {
    page = Math.max(0, page - 1);
    render();
  });

  next.addEventListener('click', () => {
    page = page + 1;
    render();
  });

  renderHead();
  render();
"""


def _build_data_table_page(frame, page_size, csv_name):
    """Build the standalone, client side filterable/sortable table page."""
    import json

    payload = frame.to_json(orient="split", index=False)

    # '<' can only appear inside JSON strings, so escaping it keeps the payload
    # from ever closing the script tag early
    payload = payload.replace("<", "\\u003c")

    sizes = sorted({25, 50, 100, 200, 500, 1000, int(page_size)})
    options = "".join(
        f'<option value="{size}"{" selected" if size == int(page_size) else ""}>'
        f"{size} rows</option>"
        for size in sizes
    )

    return (
        "<style>" + _BASE_CSS + _TABLE_CSS + "</style>"
        '<div class="eapm-wrap">'
        "<h3>Relax data</h3>"
        f"<p>Every pose of the calculation. The same table is written to "
        f"<code>{csv_name}</code>.</p>"
        '<div class="eapm-controls">'
        '<input id="eapm-search" type="search" placeholder="Filter rows..."/>'
        f'<select id="eapm-page-size">{options}</select>'
        '<span class="eapm-count" id="eapm-count"></span>'
        "</div>"
        '<div class="eapm-scroll" style="max-height:60vh">'
        '<table class="eapm-table">'
        '<thead id="eapm-head"></thead><tbody id="eapm-body"></tbody>'
        "</table></div>"
        '<div class="eapm-pager">'
        '<button id="eapm-prev" type="button">Previous</button>'
        '<button id="eapm-next" type="button">Next</button>'
        '<span class="eapm-count" id="eapm-page-info"></span>'
        "</div></div>"
        "<script>const PAYLOAD = " + payload + ";"
        "const PAGE_SIZE = " + str(int(page_size)) + ";</script>"
        "<script>" + _TABLE_JS + "</script>"
    )


def _run_extraction_jobs(jobs, cpus, verbose, interpreter):
    """Generate the parallel scripts and run them to completion."""
    import os
    import re
    import shlex
    import subprocess

    import bsc_calculations

    script_name = "relax_commands"

    # Clean up scripts from a previous run so we never execute a stale set.
    for existing in os.listdir("."):
        if re.match(rf"^{re.escape(script_name)}(_\d+)?$", existing):
            os.remove(existing)

    bsc_calculations.local.parallel(jobs, cpus=cpus, script_name=script_name)

    numbered = sorted(
        f for f in os.listdir(".") if re.match(rf"^{re.escape(script_name)}_\d+$", f)
    )
    if not numbered:
        raise ValueError("No extraction scripts were generated.")

    # The generated commands call a bare 'python', which is not necessarily on
    # the PATH of the server process and, when it is, is not necessarily the
    # environment holding PyRosetta. Point them at the configured interpreter.
    quoted = shlex.quote(interpreter)
    for script in numbered:
        with open(script) as sf:
            contents = sf.read()

        patched = re.sub(r"^python(?=\s)", quoted, contents, flags=re.MULTILINE)

        with open(script, "w") as sf:
            sf.write(patched)

    print(f"Running {len(jobs)} extraction job(s) across {len(numbered)} script(s)...")

    output = None if verbose else subprocess.DEVNULL
    processes = [
        subprocess.Popen(["bash", script], stdout=output, stderr=subprocess.PIPE)
        for script in numbered
    ]

    failures = []
    for script, process in zip(numbered, processes):
        _, stderr = process.communicate()
        if process.returncode != 0:
            message = (stderr or b"").decode("utf-8", "replace").strip()
            failures.append(f"{script} (exit {process.returncode}): {message[:300]}")

    if failures:
        raise ValueError(
            "Some Rosetta extraction scripts failed:\n  " + "\n  ".join(failures)
        )

    print("All extraction jobs finished.")


def _build_score_plot(relax_data, lengths, plot_type, max_models, output_path):
    """
    Draw the per-model score distributions, ordered best (lowest) mean score first.

    When ``lengths`` is given a second panel shows the score normalised per residue,
    which is the quantity the notebook uses for the Pareto ranking.

    Returns a note describing any truncation, or an empty string.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    data = relax_data.reset_index()
    if "Model" not in data.columns:
        raise ValueError("The Rosetta data has no 'Model' index level to group by.")

    panels = [("score", "Rosetta score (kcal/mol)")]

    if lengths:
        data["_length"] = data["Model"].map(lengths)
        missing = sorted(set(data.loc[data["_length"].isna(), "Model"]))
        if missing:
            print(f"Warning: {len(missing)} model(s) have no length and are excluded "
                  f"from the per-residue panel, e.g. {missing[:3]}.")
        with_length = data.dropna(subset=["_length"]).copy()
        if not with_length.empty:
            with_length["score_per_residue"] = (
                with_length["score"] / with_length["_length"]
            )
            panels.append(("score_per_residue", "Rosetta score per residue"))
    else:
        with_length = None

    # Order models by mean score, best (most negative) first
    order = data.groupby("Model")["score"].mean().sort_values().index.tolist()
    note = ""
    if max_models and len(order) > max_models:
        note = (f"Showing the {max_models} best-scoring models out of {len(order)}.")
        order = order[:max_models]
        print(note)

    width = max(10, min(0.28 * len(order), 40))
    fig, axes = plt.subplots(
        len(panels), 1, figsize=(width, 5.5 * len(panels)), squeeze=False
    )

    for index, (column, label) in enumerate(panels):
        axis = axes[index][0]
        frame = data if column == "score" else with_length
        subset = frame[frame["Model"].isin(order)]

        if plot_type == "box":
            sns.boxplot(data=subset, x="Model", y=column, order=order, ax=axis)
        elif plot_type == "strip":
            sns.stripplot(data=subset, x="Model", y=column, order=order, ax=axis,
                          size=3, alpha=0.6)
        else:
            sns.violinplot(data=subset, x="Model", y=column, order=order,
                           inner="quartile", ax=axis)

        axis.set_ylabel(label)
        axis.set_xlabel("")
        axis.grid(True, axis="y", alpha=0.3)
        axis.tick_params(axis="x", labelrotation=90, labelsize=6)

    axes[0][0].set_title(
        "Rosetta relax score distributions (ordered by mean score)", fontsize=14
    )
    axes[-1][0].set_xlabel("Model")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    return note


def analyse_rosetta_relax(block: PluginBlock):
    """
    Extract and summarise the Rosetta relax scores.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import json
    import os
    import sys

    import pandas as pd
    import prepare_proteins

    # pylint: enable=import-outside-toplevel

    if block.remote.name != "Local":
        raise ValueError("This block only works on the local machine.")

    models_folder = block.inputs.get(modelsFolder.id, None)
    if not models_folder or not os.path.isdir(models_folder):
        raise ValueError("A valid models folder must be provided.")

    relax_folder = block.inputs.get(relaxFolder.id, None)
    if not relax_folder or not os.path.isdir(relax_folder):
        raise ValueError("A valid Rosetta relax folder must be provided.")

    cpus = int(block.variables.get(cpusVariable.id, 4) or 4)
    skip_finished = block.variables.get(skipFinishedVariable.id, True)
    overwrite = block.variables.get(overwriteVariable.id, False)
    energy_by_residue = block.variables.get(energyByResidueVariable.id, False)
    verbose = block.variables.get(verboseVariable.id, False)
    max_rows = int(block.variables.get(maxRowsVariable.id, 200) or 200)
    plot_type = block.variables.get(plotTypeVariable.id) or "violin"
    max_models_plot = int(block.variables.get(maxModelsPlotVariable.id, 60) or 0)

    lengths_path = block.inputs.get(lengthsFile.id, None)
    lengths = {}
    if lengths_path and os.path.isfile(lengths_path):
        with open(lengths_path) as jf:
            lengths = {str(k): float(v) for k, v in json.load(jf).items()}
        print(f"Loaded {len(lengths)} sequence lengths for per-residue normalisation.")

    print("Loading models...")
    models = prepare_proteins.proteinModels(models_folder, ignore_biopython_warnings=True)

    analysis_kwargs = {
        "energy_by_residue": energy_by_residue,
        "overwrite": overwrite,
        "verbose": verbose,
    }

    # ---- Phase 1: build and run the extraction jobs ----
    jobs = models.analyseRosettaCalculation(
        relax_folder, return_jobs=True, skip_finished=skip_finished, **analysis_kwargs
    )

    if jobs:
        interpreter = block.config.get("pyrosetta_python") or sys.executable
        _run_extraction_jobs(jobs, cpus, verbose, interpreter)
    else:
        print("No extraction jobs pending; reading existing results.")

    # ---- Phase 2: collect the data ----
    print("Collecting Rosetta relax data...")
    relax_data = models.analyseRosettaCalculation(relax_folder, **analysis_kwargs)

    if relax_data is None or len(relax_data) == 0:
        raise ValueError(
            f"No Rosetta relax data was produced from '{relax_folder}'. "
            "Check that the relax calculation completed."
        )

    if "score" not in relax_data.columns:
        raise ValueError(
            "The Rosetta data has no 'score' column. Available columns: "
            f"{list(relax_data.columns)}"
        )

    relax_data_output = "rosetta_relax_data.csv"
    relax_data.to_csv(relax_data_output)

    # Per-model score statistics. 'Model' is an index level of the returned frame.
    grouped = relax_data.groupby("Model")["score"]
    summary = grouped.agg(["mean", "min", "max", "std", "count"])
    summary = summary.rename(
        columns={
            "mean": "score_mean",
            "min": "score_min",
            "max": "score_max",
            "std": "score_std",
            "count": "poses",
        }
    ).sort_values("score_mean")

    if lengths:
        summary["length"] = [lengths.get(str(model)) for model in summary.index]
        summary["score_per_residue"] = summary["score_mean"] / summary["length"]

    summary_output = "rosetta_score_summary.csv"
    summary.to_csv(summary_output)

    mean_scores = {str(model): float(value) for model, value in grouped.mean().items()}
    mean_scores_output = "rosetta_mean_scores.json"
    with open(mean_scores_output, "w") as jf:
        json.dump(mean_scores, jf, indent=2)

    print(
        f"Analysed {len(summary)} model(s) / {len(relax_data)} pose(s). "
        f"Best mean score: {summary['score_mean'].iloc[0]:.2f} "
        f"({summary.index[0]})"
    )

    # ---- Score plot ----
    import base64 as _base64

    plot_output = "rosetta_score_plot.png"
    plot_note = _build_score_plot(
        relax_data, lengths, plot_type, max_models_plot, plot_output
    )
    with open(plot_output, "rb") as img:
        plot_encoded = _base64.b64encode(img.read()).decode("ascii")
    plot_html = (
        f"<h3>Score distributions</h3>"
        f"{f'<p><em>{plot_note}</em></p>' if plot_note else ''}"
        f'<div style="overflow-x:auto; margin-bottom:24px">'
        f'<img src="data:image/png;base64,{plot_encoded}" style="max-width:100%"/></div>'
    )

    # ---- Results views ----
    # Two separate pages: the plot plus the per-model summary, and the full
    # per-pose data as its own interactive table.
    pd.set_option("display.max_colwidth", 200)

    summary_html = summary.round(3).to_html(classes="eapm-table", border=0)

    summary_page = f"""
    <style>{_BASE_CSS}</style>
    <div class="eapm-wrap">
      <h3>Rosetta relax summary</h3>
      <p>{len(summary)} models &middot; {len(relax_data)} poses &middot;
         best mean score <b>{summary['score_mean'].iloc[0]:.2f}</b>
         ({summary.index[0]})</p>
      {plot_html}
      <h3>Per-model summary</h3>
      <div class="eapm-scroll">{summary_html}</div>
    </div>
    """
    Extensions().loadHTML(summary_page, title="Rosetta relax summary")

    data_page = _build_data_table_page(
        relax_data.reset_index().round(3), max_rows, relax_data_output
    )
    Extensions().loadHTML(data_page, title="Rosetta relax data")

    block.setOutput(relaxDataFile.id, relax_data_output)
    block.setOutput(scoreSummaryFile.id, summary_output)
    block.setOutput(meanScoresFile.id, mean_scores_output)
    block.setOutput(plotFile.id, plot_output)


analyseRosettaRelaxBlock = PluginBlock(
    category="Rosetta",
    name="Analyse Rosetta Relax",
    id="analyse_rosetta_relax",
    description="Extract the Rosetta relax scores, summarise them per model and "
    "display them in a table view.",
    inputs=[modelsFolder, relaxFolder, lengthsFile],
    variables=[
        cpusVariable,
        skipFinishedVariable,
        overwriteVariable,
        energyByResidueVariable,
        verboseVariable,
        maxRowsVariable,
        plotTypeVariable,
        maxModelsPlotVariable,
    ],
    outputs=[relaxDataFile, scoreSummaryFile, meanScoresFile, plotFile],
    action=analyse_rosetta_relax,
)
