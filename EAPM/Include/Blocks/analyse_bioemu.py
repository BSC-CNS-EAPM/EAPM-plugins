"""
Module containing the BioEmu analysis block for the EAPM plugin.

Reproduces the UPO notebook's BioEmu analysis (cells 41-50) in a single block::

    rmsd = selected_sequences.computeBioEmuRMSD(bioemu_folder, models_paths)
    rmsf = selected_sequences.computeBioEmuRMSF(bioemu_folder, models_paths, plot=False)
    native_contacts = selected_sequences.computeNativeContacts('native_contacts',
                                                              bioemu_folder, models_paths)
    Q  = selected_sequences.computeFractionOfNativeContacts(native_contacts)
    Gf = selected_sequences.computeFoldingFreeEnergy(Q)

The four metrics are merged into one per-model table using the notebook's column
names (cell 55), so it can be fed straight into the Pareto ranking, and the
distributions are drawn in a single multi-panel view.

Native contacts require SMOG; that phase can be switched off, in which case only
the RMSD and RMSF metrics are produced.
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
sequencesFile = PluginVariable(
    id="sequences_file",
    name="Sequences file",
    description="FASTA (or JSON {name: sequence}) with the sampled sequences.",
    type=VariableTypes.FILE,
    allowedValues=["fasta", "fa", "faa", "json"],
)
bioemuFolder = PluginVariable(
    id="bioemu_folder",
    name="BioEmu folder",
    description="BioEmu sampling folder produced by the 'BioEmu Sampling' block.",
    type=VariableTypes.FOLDER,
)
modelsFolder = PluginVariable(
    id="models_folder",
    name="Reference models folder",
    description="Folder with the reference structures (e.g. 'selected_models'), used "
    "as the RMSD/RMSF reference and as the native state for the contacts.",
    type=VariableTypes.FOLDER,
)

# ==========================#
# Variables (parameters)
# ==========================#
computeContactsVariable = PluginVariable(
    id="compute_native_contacts",
    name="Compute native contacts",
    description="Compute the fraction of native contacts (Q) and the folding free "
    "energy. Requires SMOG; disable to get only RMSD and RMSF.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)
contactsFolderVariable = PluginVariable(
    id="native_contacts_folder",
    name="Native contacts folder",
    description="Working folder for the native-contacts calculation.",
    type=VariableTypes.STRING,
    defaultValue="native_contacts",
)
qMethodVariable = PluginVariable(
    id="q_method",
    name="Q method",
    description="Method used to decide whether a native contact is formed.",
    type=VariableTypes.STRING,
    defaultValue="hard",
)
inflationVariable = PluginVariable(
    id="inflation",
    name="Contact inflation",
    description="Distance inflation factor applied to the native contact cutoff.",
    type=VariableTypes.FLOAT,
    defaultValue=1.2,
)
relToleranceVariable = PluginVariable(
    id="rel_tolerance",
    name="Relative tolerance",
    description="Relative tolerance used by the soft contact criterion.",
    type=VariableTypes.FLOAT,
    defaultValue=0.2,
)
qThresholdVariable = PluginVariable(
    id="q_threshold",
    name="Q threshold",
    description="Q value above which a conformation counts as folded.",
    type=VariableTypes.FLOAT,
    defaultValue=0.8,
)
temperatureVariable = PluginVariable(
    id="temperature",
    name="Temperature (K)",
    description="Temperature used for the folding free energy.",
    type=VariableTypes.FLOAT,
    defaultValue=300.0,
)
plotTypeVariable = PluginVariable(
    id="plot_type",
    name="Plot type",
    description="How to show the per-model RMSD and Q distributions.",
    type=VariableTypes.STRING,
    defaultValue="violin",
    allowedValues=["violin", "box"],
)
maxModelsPlotVariable = PluginVariable(
    id="max_models_plot",
    name="Max models in plot",
    description="Only this many models are drawn, keeping the x-axis readable. "
    "Set to 0 to plot every model.",
    type=VariableTypes.INTEGER,
    defaultValue=60,
)

# ==========================#
# Variable outputs
# ==========================#
metricsFile = PluginVariable(
    id="metrics_file",
    name="BioEmu metrics",
    description="CSV with the combined per-model metrics (RMSD, RMSF, Q and dGf).",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
rmsdFile = PluginVariable(
    id="rmsd_file",
    name="RMSD values",
    description="JSON mapping each model to its per-sample RMSD values.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
rmsfFile = PluginVariable(
    id="rmsf_file",
    name="RMSF values",
    description="JSON mapping each model to its per-residue RMSF values.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
qFile = PluginVariable(
    id="q_file",
    name="Fraction of native contacts",
    description="CSV with the per-model/pose fraction of native contacts (Q).",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)
plotFile = PluginVariable(
    id="plot_file",
    name="Metrics plot",
    description="PNG with the RMSD, RMSF, Q and folding free energy panels.",
    type=VariableTypes.FILE,
    allowedValues=["png"],
)


def _ensure_fasta(sequences_path: str) -> str:
    """Return a FASTA path, converting a JSON ``{name: sequence}`` file if needed."""
    import json
    import os

    if not sequences_path.lower().endswith(".json"):
        return sequences_path

    with open(sequences_path) as jf:
        data = json.load(jf)
    if not isinstance(data, dict):
        raise ValueError("The JSON sequences file must map sequence names to sequences.")

    fasta_path = os.path.join(os.getcwd(), "bioemu_analysis_input.fasta")
    with open(fasta_path, "w") as of:
        for name, sequence in data.items():
            of.write(f">{name}\n{sequence}\n")
    return fasta_path


def _to_list(values):
    """Normalise a numpy array / Series / iterable of numbers into a plain list."""
    if values is None:
        return []
    if hasattr(values, "tolist"):
        return [float(v) for v in values.tolist()]
    return [float(v) for v in values]


def _build_metrics(rmsd, rmsf, q_frame, gf_frame):
    """Merge the BioEmu metrics into one per-model table (notebook cell 55 naming)."""
    import numpy as np
    import pandas as pd

    frames = []

    if rmsd:
        frames.append(
            pd.DataFrame(
                [
                    {
                        "Model": model,
                        "RMSD_Average": float(np.mean(values)) if len(values) else np.nan,
                        "RMSD_StdDev": float(np.std(values)) if len(values) else np.nan,
                    }
                    for model, values in rmsd.items()
                ]
            ).set_index("Model")
        )

    if rmsf:
        frames.append(
            pd.DataFrame(
                [
                    {
                        "Model": model,
                        "RMSF_Average": float(np.mean(values)) if len(values) else np.nan,
                    }
                    for model, values in rmsf.items()
                ]
            ).set_index("Model")
        )

    if q_frame is not None and len(q_frame):
        q_stats = q_frame.groupby("Model")["Q"].agg(["mean", "std"])
        q_stats.columns = ["Q_Average", "Q_StdDev"]
        frames.append(q_stats)

    if gf_frame is not None and len(gf_frame):
        column = next(
            (c for c in gf_frame.columns if "G" in c and "f" in c), gf_frame.columns[0]
        )
        frames.append(gf_frame[[column]].rename(columns={column: "Folding_FreeEnergy_dGf"}))

    if not frames:
        raise ValueError("No BioEmu metrics could be computed.")

    return pd.concat(frames, axis=1)


def _build_plot(rmsd, rmsf, q_frame, gf_frame, plot_type, max_models, output_path):
    """Draw the RMSD / RMSF / Q / dGf panels, returning a truncation note."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns

    # Order models by mean RMSD (most native-like first)
    if rmsd:
        order = sorted(rmsd, key=lambda m: float(np.mean(rmsd[m])) if len(rmsd[m]) else np.inf)
    elif rmsf:
        order = sorted(rmsf, key=lambda m: float(np.mean(rmsf[m])))
    else:
        order = sorted(q_frame.index.get_level_values("Model").unique())

    note = ""
    if max_models and len(order) > max_models:
        note = f"Showing {max_models} of {len(order)} models."
        order = order[:max_models]

    panels = []
    if rmsd:
        panels.append("rmsd")
    if rmsf:
        panels.append("rmsf")
    if q_frame is not None and len(q_frame):
        panels.append("q")
    if gf_frame is not None and len(gf_frame):
        panels.append("gf")

    width = max(10, min(0.28 * len(order), 40))
    fig, axes = plt.subplots(len(panels), 1, figsize=(width, 4.8 * len(panels)), squeeze=False)

    def _dist(axis, frame, column, label):
        if plot_type == "box":
            sns.boxplot(data=frame, x="Model", y=column, order=order, ax=axis)
        else:
            sns.violinplot(data=frame, x="Model", y=column, order=order,
                           inner="quartile", ax=axis)
        axis.set_ylabel(label)

    for index, panel in enumerate(panels):
        axis = axes[index][0]

        if panel == "rmsd":
            frame = pd.DataFrame(
                [{"Model": m, "RMSD": v} for m in order for v in _to_list(rmsd.get(m, []))]
            )
            _dist(axis, frame, "RMSD", "RMSD (Å)")
        elif panel == "rmsf":
            frame = pd.DataFrame(
                [
                    {"Model": m, "RMSF": float(np.mean(rmsf[m]))}
                    for m in order
                    if m in rmsf
                ]
            )
            sns.barplot(data=frame, x="Model", y="RMSF", order=order, ax=axis)
            axis.set_ylabel("Average RMSF (Å)")
        elif panel == "q":
            frame = q_frame.reset_index()
            frame = frame[frame["Model"].isin(order)]
            _dist(axis, frame, "Q", "Fraction of native contacts (Q)")
        else:
            column = next(
                (c for c in gf_frame.columns if "G" in c and "f" in c), gf_frame.columns[0]
            )
            frame = gf_frame.reset_index()
            frame.columns = ["Model"] + list(frame.columns[1:])
            frame = frame[frame["Model"].isin(order)]
            sns.barplot(data=frame, x="Model", y=column, order=order, ax=axis)
            axis.axhline(0, color="black", linewidth=0.8, linestyle="--")
            axis.set_ylabel("ΔG_f (kcal/mol)")

        axis.set_xlabel("")
        axis.grid(True, axis="y", alpha=0.3)
        axis.tick_params(axis="x", labelrotation=90, labelsize=6)

    axes[0][0].set_title("BioEmu ensemble metrics (ordered by mean RMSD)", fontsize=14)
    axes[-1][0].set_xlabel("Model")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)

    return note


def analyse_bioemu(block: PluginBlock):
    """
    Compute the BioEmu ensemble metrics and display them.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import base64
    import json
    import os

    import prepare_proteins

    # pylint: enable=import-outside-toplevel

    if block.remote.name != "Local":
        raise ValueError("This block only works on the local machine.")

    sequences_path = block.inputs.get(sequencesFile.id, None)
    if not sequences_path or not os.path.isfile(sequences_path):
        raise ValueError("A valid sequences file (FASTA or JSON) must be provided.")

    bioemu_folder = block.inputs.get(bioemuFolder.id, None)
    if not bioemu_folder or not os.path.isdir(bioemu_folder):
        raise ValueError("A valid BioEmu folder must be provided.")

    models_folder = block.inputs.get(modelsFolder.id, None)
    if not models_folder or not os.path.isdir(models_folder):
        raise ValueError("A valid reference models folder must be provided.")

    compute_contacts = block.variables.get(computeContactsVariable.id, True)
    contacts_folder = block.variables.get(contactsFolderVariable.id) or "native_contacts"
    q_method = block.variables.get(qMethodVariable.id) or "hard"
    inflation = block.variables.get(inflationVariable.id, 1.2)
    rel_tolerance = block.variables.get(relToleranceVariable.id, 0.2)
    q_threshold = block.variables.get(qThresholdVariable.id, 0.8)
    temperature = block.variables.get(temperatureVariable.id, 300.0)
    plot_type = block.variables.get(plotTypeVariable.id) or "violin"
    max_models_plot = int(block.variables.get(maxModelsPlotVariable.id, 60) or 0)

    sequences_path = _ensure_fasta(sequences_path)

    print("Loading sequences and reference models...")
    sequences = prepare_proteins.sequenceModels(sequences_path)
    models = prepare_proteins.proteinModels(models_folder, ignore_biopython_warnings=True)
    models_paths = models.models_paths

    # ---- RMSD ----
    print("Computing BioEmu RMSD...")
    rmsd_raw = sequences.computeBioEmuRMSD(bioemu_folder, models_paths, plot=False)
    rmsd = {str(model): _to_list(values) for model, values in dict(rmsd_raw).items()}
    print(f"  RMSD computed for {len(rmsd)} model(s).")

    # ---- RMSF ----
    print("Computing BioEmu RMSF...")
    rmsf_raw = sequences.computeBioEmuRMSF(
        bioemu_folder, models_paths, plot=False, plot_legend=False
    )
    rmsf = {str(model): _to_list(values) for model, values in dict(rmsf_raw).items()}
    print(f"  RMSF computed for {len(rmsf)} model(s).")

    # ---- Native contacts -> Q -> dGf ----
    q_frame = None
    gf_frame = None
    if compute_contacts:
        print("Computing native contacts (requires SMOG)...")
        try:
            native_contacts = sequences.computeNativeContacts(
                contacts_folder, bioemu_folder, models_paths
            )
            q_frame = sequences.computeFractionOfNativeContacts(
                native_contacts,
                method=q_method,
                inflation=inflation,
                rel_tolerance=rel_tolerance,
            )
            gf_frame = sequences.computeFoldingFreeEnergy(
                q_frame, q_threshold=q_threshold, temperature=temperature
            )
            print(f"  Q computed for {len(q_frame)} model/pose rows.")
        except Exception as error:  # pylint: disable=broad-except
            raise ValueError(
                f"The native-contacts step failed: {error}. This phase needs SMOG; "
                "disable 'Compute native contacts' to get only RMSD and RMSF."
            ) from error
    else:
        print("Skipping native contacts (disabled); Q and dGf will not be computed.")

    # ---- Combined metrics ----
    metrics = _build_metrics(rmsd, rmsf, q_frame, gf_frame)
    metrics_output = "bioemu_metrics.csv"
    metrics.to_csv(metrics_output)
    print(f"Combined metrics for {len(metrics)} model(s): {list(metrics.columns)}")

    rmsd_output = "bioemu_rmsd.json"
    with open(rmsd_output, "w") as jf:
        json.dump(rmsd, jf, indent=2)

    rmsf_output = "bioemu_rmsf.json"
    with open(rmsf_output, "w") as jf:
        json.dump(rmsf, jf, indent=2)

    q_output = None
    if q_frame is not None:
        q_output = "bioemu_q.csv"
        q_frame.to_csv(q_output)

    # ---- Plot + view ----
    plot_output = "bioemu_metrics_plot.png"
    note = _build_plot(rmsd, rmsf, q_frame, gf_frame, plot_type, max_models_plot, plot_output)

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
      <h3>BioEmu ensemble metrics</h3>
      <p>{len(metrics)} models &middot; metrics: {", ".join(metrics.columns)}</p>
      {f'<p><em>{note}</em></p>' if note else ''}
      <div style="overflow-x:auto; margin-bottom:24px">
        <img src="data:image/png;base64,{encoded}" style="max-width:100%"/>
      </div>
      <h3>Per-model metrics</h3>
      <div class="eapm-scroll">{metrics.round(3).to_html(classes="eapm-table", border=0)}</div>
    </div>
    """
    Extensions().loadHTML(html, title="BioEmu metrics")

    block.setOutput(metricsFile.id, metrics_output)
    block.setOutput(rmsdFile.id, rmsd_output)
    block.setOutput(rmsfFile.id, rmsf_output)
    block.setOutput(plotFile.id, plot_output)
    if q_output:
        block.setOutput(qFile.id, q_output)


analyseBioEmuBlock = PluginBlock(
    category="BioEmu",
    name="Analyse BioEmu",
    id="analyse_bioemu",
    description="Compute the RMSD, RMSF, fraction of native contacts (Q) and folding "
    "free energy of a BioEmu ensemble, and display them in one table and plot.",
    inputs=[sequencesFile, bioemuFolder, modelsFolder],
    variables=[
        computeContactsVariable,
        contactsFolderVariable,
        qMethodVariable,
        inflationVariable,
        relToleranceVariable,
        qThresholdVariable,
        temperatureVariable,
        plotTypeVariable,
        maxModelsPlotVariable,
    ],
    outputs=[metricsFile, rmsdFile, rmsfFile, qFile, plotFile],
    action=analyse_bioemu,
)
