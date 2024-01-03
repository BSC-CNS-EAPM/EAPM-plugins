import os
import shutil

from HorusAPI import PluginBlock, PluginVariable, VariableTypes, VariableGroup

resultsFolderAF = PluginVariable(
    name="Models folder",
    id="results_folder",
    description="Folder containing the Alphafold output PDBs",
    type=VariableTypes.FOLDER,
    defaultValue="alphafold",
)

alphafoldOutputResultsVariable = PluginVariable(
    id="alphafold_output_results",
    name="Alphafold results",
    description="BSC JOB Alphafold results",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_results"],
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
    if block.selectedInputGroup == "alphafold_folder_group":
        models_folder = block.inputs.get("results_folder", None)
    elif block.selectedInputGroup == "alphafold_job_group":
        models_folder = block.inputs.get("alphafold_output_results", None)

        models_folder = models_folder.get("output_models", None)
    else:
        raise Exception("No valid input group selected")

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
    inputGroups=[
        VariableGroup(
            id="alphafold_folder_group",
            name="Alphafold folder",
            description="Alphafold models folder",
            variables=[
                resultsFolderAF,
            ],
        ),
        VariableGroup(
            id="alphafold_job_group",
            name="Alphafold output",
            description="Output models from BSC job calculation from Alphafold",
            variables=[
                alphafoldOutputResultsVariable,
            ],
        ),
    ],
    variables=[
        confidenceThresholdAF,
    ],
    outputs=[
        outputPDBAF,
        trimmedModelsOutputAF,
    ],
    action=trimAlphaFoldModels,
)
