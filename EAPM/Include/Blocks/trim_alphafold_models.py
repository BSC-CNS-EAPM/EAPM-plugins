"""
Module containing the trim alphafold models block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

resultsFolderAF = PluginVariable(
    name="Models",
    id="results_folder",
    description="Folder containing the Alphafold output PDBs",
    type=VariableTypes.FOLDER,
    defaultValue="alphafold",
)

confidenceThresholdAF = PluginVariable(
    name="Confidence threshold",
    id="confidence_threshold",
    description="Threshold confidence indicates the maximum"
    "confidence score at which to stop the trimming of terminal regions.",
    type=VariableTypes.FLOAT,
    defaultValue=90.0,
)

outputPDBAF = PluginVariable(
    name="First PDB",
    id="out_pdb",
    description="First PDB of the Alphafold.",
    type=VariableTypes.FILE,
    allowedValues=["pdb"],
)

trimmedModelsOutputAF = PluginVariable(
    name="Trimmed models",
    id="trimmed_models",
    description="Folder containing the trimmed models",
    type=VariableTypes.FOLDER,
    defaultValue="trimmed_models",
)


def trim_alphafold_models(block: PluginBlock):
    """
    Trims AlphaFold models based on confidence scores.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.

    Raises:
        ValueError: If no models folder is selected.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os
    import shutil

    import prepare_proteins

    # pylint: enable=import-outside-toplevel
    # Get the models folder
    models_folder = block.inputs.get(resultsFolderAF.id, None)

    if models_folder is None:
        raise ValueError("No models folder selected")

    # Create the untrimmed_models folder
    untrimmed_folder = os.path.join(os.getcwd(), "untrimmed_models")
    if not os.path.isdir(untrimmed_folder):
        os.mkdir(untrimmed_folder)

    # Copy each alphafold output model (from a specific rank) into the models folder
    rank = 0
    for model in os.listdir(models_folder):
        ranked_file = os.path.join(models_folder, model, "ranked_" + str(rank) + ".pdb")
        new_model_filename = model + ".pdb"
        if os.path.exists(ranked_file):
            shutil.copyfile(
                ranked_file,
                os.path.join(untrimmed_folder, new_model_filename),
            )

    print(f"Loading models from {untrimmed_folder}...")
    models = prepare_proteins.proteinModels(untrimmed_folder)

    print("Trimming the models...")
    confidence_threshold = float(block.variables.get(confidenceThresholdAF.id, 90))

    models.removeTerminiByConfidenceScore(confidence_threshold)

    print("Saving trimmed models")
    models.saveModels("trimmed_models")

    # Remove the untrimmed models folder
    shutil.rmtree(untrimmed_folder)

    print("Setting output of block to the results directory...")

    trimmed_folder = os.path.join(os.getcwd(), "trimmed_models")

    # Set the output
    block.setOutput(trimmedModelsOutputAF.id, trimmed_folder)

    out_pdb = None
    for file in os.listdir(trimmed_folder):
        if file.endswith(".pdb"):
            out_pdb = os.path.join(trimmed_folder, file)
            break

    block.setOutput(outputPDBAF.id, out_pdb)


trimAlphaFoldModelsBlock = PluginBlock(
    name="Trim Alphafold models",
    id="trim_alphafold_models",
    description="Trim the Alphafold models",
    inputs=[resultsFolderAF],
    variables=[
        confidenceThresholdAF,
    ],
    outputs=[
        outputPDBAF,
        trimmedModelsOutputAF,
    ],
    action=trim_alphafold_models,
)
