"""
Module containing the selected-PDB collection block for the EAPM plugin.

Reproduces the UPO notebook's "Load models" step (cells 24-25), which gathers the
structures of the selected models into a dedicated working folder::

    models_folder = 'selected_models'
    os.makedirs(models_folder, exist_ok=True)
    for model in selected_sequences:
        model_pdb = models_folder+'/'+model+'.pdb'
        if not os.path.exists(model_pdb):
            shutil.copyfile('AF_truncated/'+model+'.pdb', model_pdb)

The resulting folder is what downstream structure-based steps (Rosetta relax,
BioEmu RMSD/RMSF, ...) operate on.
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
sourceFolder = PluginVariable(
    id="source_folder",
    name="Source structures folder",
    description="Folder holding the structures of all models (e.g. 'AF_truncated').",
    type=VariableTypes.FOLDER,
)
modelsFile = PluginVariable(
    id="models_file",
    name="Selected models",
    description="The selection to collect. A JSON list of model names, a JSON "
    "{name: sequence} mapping, or a FASTA file.",
    type=VariableTypes.FILE,
    allowedValues=["json", "fasta", "fa", "faa"],
)

# ==========================#
# Variables (parameters)
# ==========================#
outputFolderVariable = PluginVariable(
    id="output_folder",
    name="Output folder",
    description="Name of the folder the selected structures are copied into.",
    type=VariableTypes.STRING,
    defaultValue="selected_models",
)
extensionVariable = PluginVariable(
    id="extension",
    name="Structure extension",
    description="Extension of the structure files in the source folder.",
    type=VariableTypes.STRING,
    defaultValue="pdb",
    allowedValues=["pdb", "cif", "mae"],
)
overwriteVariable = PluginVariable(
    id="overwrite",
    name="Overwrite existing",
    description="Re-copy structures that are already present in the output folder. "
    "Disabled by default, matching the notebook's skip-if-exists behaviour.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
ignoreMissingVariable = PluginVariable(
    id="ignore_missing",
    name="Ignore missing structures",
    description="Warn and continue when a selected model has no structure in the "
    "source folder. Disabled by default, so a missing structure is an error.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)

# ==========================#
# Variable outputs
# ==========================#
modelsFolderOutput = PluginVariable(
    id="models_folder",
    name="Selected models folder",
    description="Folder containing the collected structures of the selected models.",
    type=VariableTypes.FOLDER,
)
collectedModelsFile = PluginVariable(
    id="collected_models_file",
    name="Collected models",
    description="JSON list of the models present in the output folder.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
missingModelsFile = PluginVariable(
    id="missing_models_file",
    name="Missing models",
    description="JSON list of selected models that had no structure in the source "
    "folder (only written when 'Ignore missing structures' is enabled).",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)


def _read_model_names(models_path: str) -> list:
    """Extract the ordered list of model names from a JSON list/dict or a FASTA file."""
    import json

    if models_path.lower().endswith(".json"):
        with open(models_path) as jf:
            data = json.load(jf)
        if isinstance(data, list):
            return [str(name) for name in data]
        if isinstance(data, dict):
            return [str(name) for name in data]
        raise ValueError("The JSON models file must be a list of names or a "
                         "{name: sequence} mapping.")

    names = []
    with open(models_path) as ff:
        for line in ff:
            if line.startswith(">"):
                header = line[1:].strip().split()
                if header:
                    names.append(header[0])
    if not names:
        raise ValueError("No sequence names found in the FASTA models file.")
    return names


def collect_selected_pdbs(block: PluginBlock):
    """
    Copy the structures of the selected models into a dedicated folder.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import json
    import os
    import shutil

    # pylint: enable=import-outside-toplevel

    source_folder = block.inputs.get(sourceFolder.id, None)
    if not source_folder or not os.path.isdir(source_folder):
        raise ValueError("A valid source structures folder must be provided.")

    models_path = block.inputs.get(modelsFile.id, None)
    if not models_path or not os.path.isfile(models_path):
        raise ValueError("A valid selected-models file (JSON or FASTA) must be provided.")

    output_folder = block.variables.get(outputFolderVariable.id, "selected_models")
    extension = (block.variables.get(extensionVariable.id, "pdb") or "pdb").lstrip(".")
    overwrite = block.variables.get(overwriteVariable.id, False)
    ignore_missing = block.variables.get(ignoreMissingVariable.id, False)

    model_names = _read_model_names(models_path)
    # Preserve order but drop duplicates
    seen = set()
    models = [m for m in model_names if not (m in seen or seen.add(m))]

    os.makedirs(output_folder, exist_ok=True)

    collected: list = []
    missing: list = []
    copied = 0
    skipped = 0

    for model in models:
        destination = os.path.join(output_folder, f"{model}.{extension}")
        source = os.path.join(source_folder, f"{model}.{extension}")

        if os.path.exists(destination) and not overwrite:
            collected.append(model)
            skipped += 1
            continue

        if not os.path.isfile(source):
            missing.append(model)
            continue

        shutil.copyfile(source, destination)
        collected.append(model)
        copied += 1

    if missing:
        message = (
            f"{len(missing)} selected model(s) have no '{extension}' structure in "
            f"'{source_folder}' (e.g. {missing[:3]})."
        )
        if not ignore_missing:
            raise ValueError(
                message + " Enable 'Ignore missing structures' to skip them instead."
            )
        print(f"Warning: {message}")

    print(f"Collected {len(collected)} model(s) into '{output_folder}' "
          f"({copied} copied, {skipped} already present).")

    collected_output = "collected_models.json"
    with open(collected_output, "w") as jf:
        json.dump(collected, jf, indent=2)

    block.setOutput(modelsFolderOutput.id, output_folder)
    block.setOutput(collectedModelsFile.id, collected_output)

    if missing:
        missing_output = "missing_models.json"
        with open(missing_output, "w") as jf:
            json.dump(missing, jf, indent=2)
        block.setOutput(missingModelsFile.id, missing_output)


collectSelectedPDBsBlock = PluginBlock(
    category="Structure Preparation",
    name="Collect Selected Structures",
    id="collect_selected_pdbs",
    description="Copy the structures of a selected set of models from a source "
    "folder into a dedicated working folder.",
    inputs=[sourceFolder, modelsFile],
    variables=[
        outputFolderVariable,
        extensionVariable,
        overwriteVariable,
        ignoreMissingVariable,
    ],
    outputs=[modelsFolderOutput, collectedModelsFile, missingModelsFile],
    action=collect_selected_pdbs,
)
