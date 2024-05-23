import os
import shutil

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
fasta_fileAF = PluginVariable(
    name="Fasta file",
    id="fasta_file",
    description="The input fasta file.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# ==========================#
# Variables
# ==========================#
outputAF = PluginVariable(
    name="Alphafold simulation folder",
    id="folder_name",
    description="The name of the folder where the simulation will be stored.",
    type=VariableTypes.STRING,
    defaultValue="alphafold",
)
removeExistingResults = PluginVariable(
    name="Remove existing results",
    id="remove_existing_results",
    description="Remove existing results",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)

# Output variables
outputModelsVariable = PluginVariable(
    id="models",
    name="Alphafold models",
    description="The output models",
    type=VariableTypes.FOLDER,
)


def finalAlhafoldAction(block: PluginBlock):

    resultsFolder = "alphafold"
    downloaded_path = "/home/perry/data/acanella/testHorus/all_test/alphafold"

    output_models_folder = os.path.join(downloaded_path, resultsFolder, "output_models")

    if not os.path.exists("structures"):
        os.makedirs("structures")

    rank = 0
    for model in os.listdir(output_models_folder):
        if os.path.exists(f"{output_models_folder}/" + model + "/ranked_" + str(rank) + ".pdb"):
            shutil.copyfile(
                f"{output_models_folder}/" + model + "/ranked_" + str(rank) + ".pdb",
                "structures/" + model + ".pdb",
            )

    block.setOutput(outputModelsVariable.id, "structures")


from utils import BSC_JOB_VARIABLES

testBlock = PluginBlock(
    name="Test block",
    description="Test",
    action=finalAlhafoldAction,
    variables=BSC_JOB_VARIABLES + [outputAF, removeExistingResults],
    inputs=[fasta_fileAF],
    outputs=[outputModelsVariable],
)
