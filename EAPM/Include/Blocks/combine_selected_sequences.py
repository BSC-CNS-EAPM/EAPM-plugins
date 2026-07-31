"""
Module containing the selected-sequences combination block for the EAPM plugin.

Reproduces the UPO notebook's "Combine selected models" step (cell 23), which
merges the sets selected from the ProteinMPNN *vanilla* and *soluble* scores::

    selected_common_sequences = {**selected_vanilla_sequences, **selected_soluble_sequences}

Here the user chooses which ProteinMPNN model to take the selection from:
``vanilla``, ``soluble`` or ``both`` (the union, as in the notebook).

The block also emits the per-model *origin* ("vanilla only", "soluble only" or
"common"), which the notebook recomputes over and over in its downstream plots.
Origin is always evaluated against both input sets when they are available, so it
stays informative even when only one set is selected.
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
vanillaSequencesFile = PluginVariable(
    id="vanilla_sequences_file",
    name="Vanilla selected sequences",
    description="Sequences selected using the ProteinMPNN vanilla scores. "
    "FASTA or JSON ({name: sequence}).",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta", "fa", "faa", "json"],
)
solubleSequencesFile = PluginVariable(
    id="soluble_sequences_file",
    name="Soluble selected sequences",
    description="Sequences selected using the ProteinMPNN soluble scores. "
    "FASTA or JSON ({name: sequence}).",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta", "fa", "faa", "json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
modelSelectionVariable = PluginVariable(
    id="model_selection",
    name="ProteinMPNN model",
    description="Which selection to use: 'vanilla' or 'soluble' for a single set, "
    "or 'both' for the union of the two (the notebook's behaviour).",
    type=VariableTypes.STRING,
    defaultValue="both",
    allowedValues=["vanilla", "soluble", "both"],
)

# ==========================#
# Variable outputs
# ==========================#
combinedSequencesFile = PluginVariable(
    id="combined_sequences_file",
    name="Combined sequences (FASTA)",
    description="FASTA with the sequences of the combined selection.",
    type=VariableTypes.FILE,
    allowedValues=["fasta"],
)
combinedSequencesJson = PluginVariable(
    id="combined_sequences_json",
    name="Combined sequences (JSON)",
    description="JSON mapping each selected model to its sequence.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
selectedModelsFile = PluginVariable(
    id="selected_models_file",
    name="Selected models",
    description="JSON list with the names of the combined selection.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
originFile = PluginVariable(
    id="origin_file",
    name="Model origin",
    description="JSON mapping each selected model to its origin: "
    "'vanilla only', 'soluble only' or 'common'.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
originTableFile = PluginVariable(
    id="origin_table_file",
    name="Origin table",
    description="CSV with model, origin, and whether it appears in the vanilla "
    "and/or soluble selection.",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)


def _load_sequences(sequences_path: str) -> dict:
    """Read sequences from a FASTA or JSON file into a ``{name: sequence}`` dict."""
    import json

    if sequences_path.lower().endswith(".json"):
        with open(sequences_path) as jf:
            data = json.load(jf)
        if isinstance(data, list):
            # A plain list of model names (e.g. 'selected_models.json'): no sequences.
            return {str(name): "" for name in data}
        if not isinstance(data, dict):
            raise ValueError("The JSON sequences file must map sequence names to sequences.")
        return {str(name): str(seq) for name, seq in data.items()}

    sequences: dict = {}
    name = None
    chunks: list = []
    with open(sequences_path) as ff:
        for line in ff:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def combine_selected_sequences(block: PluginBlock):
    """
    Combine the vanilla and/or soluble selections according to the chosen model.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import csv
    import json
    import os

    # pylint: enable=import-outside-toplevel

    vanilla_path = block.inputs.get(vanillaSequencesFile.id, None)
    soluble_path = block.inputs.get(solubleSequencesFile.id, None)

    model_selection = (block.variables.get(modelSelectionVariable.id, "both") or "both").lower()
    if model_selection not in ("vanilla", "soluble", "both"):
        raise ValueError(
            f"Unknown model selection '{model_selection}'. "
            "Choose one of: vanilla, soluble, both."
        )

    has_vanilla = bool(vanilla_path) and os.path.isfile(vanilla_path)
    has_soluble = bool(soluble_path) and os.path.isfile(soluble_path)

    if model_selection == "vanilla" and not has_vanilla:
        raise ValueError("Model selection is 'vanilla' but no vanilla sequences file was given.")
    if model_selection == "soluble" and not has_soluble:
        raise ValueError("Model selection is 'soluble' but no soluble sequences file was given.")
    if model_selection == "both" and not (has_vanilla or has_soluble):
        raise ValueError(
            "Model selection is 'both' but neither a vanilla nor a soluble "
            "sequences file was given."
        )

    vanilla_sequences = _load_sequences(vanilla_path) if has_vanilla else {}
    soluble_sequences = _load_sequences(soluble_path) if has_soluble else {}

    if model_selection == "both" and not (has_vanilla and has_soluble):
        missing = "soluble" if has_vanilla else "vanilla"
        print(f"Warning: model selection is 'both' but the {missing} file is missing; "
              f"continuing with the available set only.")

    vanilla_set = set(vanilla_sequences)
    soluble_set = set(soluble_sequences)

    # Build the combined selection. For 'both' this is the notebook's
    # {**vanilla, **soluble}: vanilla entries first, then the soluble-only ones.
    if model_selection == "vanilla":
        combined = dict(vanilla_sequences)
    elif model_selection == "soluble":
        combined = dict(soluble_sequences)
    else:
        combined = {**vanilla_sequences, **soluble_sequences}

    def get_origin(model: str) -> str:
        in_vanilla = model in vanilla_set
        in_soluble = model in soluble_set
        if in_vanilla and in_soluble:
            return "common"
        if in_vanilla:
            return "vanilla only"
        if in_soluble:
            return "soluble only"
        return "unknown"

    origins = {model: get_origin(model) for model in combined}

    counts: dict = {}
    for origin in origins.values():
        counts[origin] = counts.get(origin, 0) + 1

    print(f"Model selection: {model_selection}")
    print(f"  vanilla selection: {len(vanilla_set)} models")
    print(f"  soluble selection: {len(soluble_set)} models")
    print(f"  combined: {len(combined)} models "
          f"({', '.join(f'{k}: {v}' for k, v in sorted(counts.items()))})")

    missing_sequences = [m for m, s in combined.items() if not s]
    if missing_sequences:
        print(f"Warning: {len(missing_sequences)} models have no sequence "
              f"(inputs may be name-only lists), e.g. {missing_sequences[:3]}.")

    # ---- Write outputs into the flow working directory ----
    combined_fasta = "combined_selected_sequences.fasta"
    with open(combined_fasta, "w") as of:
        for model, sequence in combined.items():
            if sequence:
                of.write(f">{model}\n{sequence}\n")
    block.setOutput(combinedSequencesFile.id, combined_fasta)

    combined_json = "combined_selected_sequences.json"
    with open(combined_json, "w") as jf:
        json.dump(combined, jf, indent=2)
    block.setOutput(combinedSequencesJson.id, combined_json)

    models_json = "combined_selected_models.json"
    with open(models_json, "w") as jf:
        json.dump(list(combined.keys()), jf, indent=2)
    block.setOutput(selectedModelsFile.id, models_json)

    origin_json = "combined_model_origin.json"
    with open(origin_json, "w") as jf:
        json.dump(origins, jf, indent=2)
    block.setOutput(originFile.id, origin_json)

    origin_csv = "combined_model_origin.csv"
    with open(origin_csv, "w", newline="") as cf:
        writer = csv.DictWriter(
            cf, fieldnames=["model", "origin", "in_vanilla", "in_soluble"]
        )
        writer.writeheader()
        for model in combined:
            writer.writerow(
                {
                    "model": model,
                    "origin": origins[model],
                    "in_vanilla": model in vanilla_set,
                    "in_soluble": model in soluble_set,
                }
            )
    block.setOutput(originTableFile.id, origin_csv)


combineSelectedSequencesBlock = PluginBlock(
    category="Clustering & Selection",
    name="Combine Selected Sequences",
    id="combine_selected_sequences",
    description="Combine the ProteinMPNN vanilla and/or soluble selections "
    "(vanilla, soluble or both) and tag each model's origin.",
    inputs=[vanillaSequencesFile, solubleSequencesFile],
    variables=[modelSelectionVariable],
    outputs=[
        combinedSequencesFile,
        combinedSequencesJson,
        selectedModelsFile,
        originFile,
        originTableFile,
    ],
    action=combine_selected_sequences,
)
