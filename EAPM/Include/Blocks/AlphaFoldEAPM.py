"""
Module containing the AlphaFold block for the EAPM plugin
"""

import os
from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

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


# Alphafold action block
def initialAlphafold(block: SlurmBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    # Loading plugin variables
    fastaFile = block.inputs["fasta_file"]
    if fastaFile == "None":
        raise Exception("No fasta file provided.")

    folderName = block.variables.get("folder_name", "alphafold")
    removeExisting = block.variables.get("remove_existing_results", False)

    # If folder already exists, raise exception
    if removeExisting and os.path.exists(folderName):
        os.system("rm -rf " + folderName)

    if not removeExisting and os.path.exists(folderName):
        raise Exception(
            "The folder {} already exists. Please, choose another name or remove it.".format(
                folderName
            )
        )

    block.extraData["folder_name"] = folderName

    import prepare_proteins

    print("Loading fasta files...")

    sequences = prepare_proteins.sequenceModels(fastaFile)

    print("Setting up AlphaFold...")

    jobs = sequences.setUpAlphaFold(folderName)

    from utils import launchCalculationAction

    launchCalculationAction(block, jobs, folderName)


def finalAlhafoldAction(block: SlurmBlock):
    from utils import downloadResultsAction

    downloaded_path = downloadResultsAction(block)

    resultsFolder = block.extraData["folder_name"]

    output_models_folder = os.path.join(downloaded_path, resultsFolder, "output_models")

    block.setOutput(outputModelsVariable.id, output_models_folder)


from utils import BSC_JOB_VARIABLES

alphafoldBlock = SlurmBlock(
    name="Alphafold",
    description="Run Alphafold. (For cte_power, marenostrum and minotauro clusters or local)",
    initialAction=initialAlphafold,
    finalAction=finalAlhafoldAction,
    variables=BSC_JOB_VARIABLES + [outputAF, removeExistingResults],
    inputs=[fasta_fileAF],
    outputs=[outputModelsVariable],
)
