import os
import shutil

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
    description="Threshold confidence indicates the maximum confidence score at which to stop the trimming of terminal regions.",
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


def trimAlphaFoldModels(block: PluginBlock):

    # Get the models folder
    models_folder = block.inputs.get("results_folder", None)

    if models_folder is None:
        raise Exception("No models folder selected")

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

    import prepare_proteins

    print(f"Loading models from {untrimmed_folder}...")
    models = prepare_proteins.proteinModels(untrimmed_folder)

    print("Trimming the models...")
    confidenceThreshold = float(block.variables.get("confidence_threshold", 90))

    models.removeTerminiByConfidenceScore(confidenceThreshold)

    print("Saving trimmed models")
    models.saveModels("trimmed_models")

    # Remove the untrimmed models folder
    shutil.rmtree(untrimmed_folder)

    print("Setting output of block to the results directory...")

    trimmed_folder = os.path.join(os.getcwd(), "trimmed_models")

    # Set the output
    block.setOutput("trimmed_models", trimmed_folder)

    outPdb = None
    for file in os.listdir(trimmed_folder):
        if file.endswith(".pdb"):
            outPdb = os.path.join(trimmed_folder, file)
            break

    block.setOutput("out_pdb", outPdb)


trimAlphaFoldModelsBlock = PluginBlock(
    name="Trim Alphafold models",
    description="Trim the Alphafold models",
    inputs=[resultsFolderAF],
    variables=[
        confidenceThresholdAF,
    ],
    outputs=[
        outputPDBAF,
        trimmedModelsOutputAF,
    ],
    action=trimAlphaFoldModels,
)
