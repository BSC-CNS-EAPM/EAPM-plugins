"""
Module containing the AlphaFold block for the EAPM plugin
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes
from utils import BSC_JOB_VARIABLES

# ==========================#
# Variable inputs
# ==========================#
fastaFile = PluginVariable(
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
output = PluginVariable(
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
def initial_alphafold(block: SlurmBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    # pylint: disable=import-outside-toplevel
    import os

    import prepare_proteins
    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel
    # Loading plugin variables
    fasta_file = block.inputs.get(fastaFile.id, None)
    if fasta_file == "None":
        raise ValueError("No fasta file provided.")

    folder_name = block.variables.get(output.id, "alphafold")
    remove_existing = block.variables.get(removeExistingResults.id, False)

    cpus_per_task = block.variables.get("cpus_per_task")
    if cpus_per_task == 1:
        print("Alphafold requires at least 20 cpus per task. Changing to 20 cpus per task.")
        block.variables["cpus_per_task"] = 20

    partiton = block.variables.get("partition")
    if partiton is None:
        print("Alphafold requires an accelerated partition. Changing to acc_bscls.")
        block.variables["partition"] = "acc_bscls"

    # If folder already exists, raise exception
    if remove_existing and os.path.exists(folder_name):
        os.system("rm -rf " + folder_name)

    if not remove_existing and os.path.exists(folder_name):
        raise ValueError(
            f"The folder {folder_name} already exists. "
            "Please, choose another name or remove it with the remove existing folder option."
        )

    block.extraData["folder_name"] = folder_name

    print("Loading fasta files...")

    sequences = prepare_proteins.sequenceModels(fasta_file)

    print("Setting up AlphaFold...")

    jobs = sequences.setUpAlphaFold(folder_name, gpu_relax=False)

    launchCalculationAction(block, jobs, "alphafold", [folder_name])


def final_alphafold(block: SlurmBlock):
    """
    Perform the final steps of the AlphaFold algorithm.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the current block.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import downloadResultsAction

    # pylint: enable=import-outside-toplevel

    downloaded_path = downloadResultsAction(block)

    results_folder = block.extraData["folder_name"]

    output_models_folder = os.path.join(downloaded_path, results_folder, "output_models")

    block.setOutput(outputModelsVariable.id, output_models_folder)


alphafoldBlock = SlurmBlock(
    name="Alphafold",
    id="Alphafold",
    description="Run Alphafold. (For marenostrum, nord3 clusters or local)",
    initialAction=initial_alphafold,
    finalAction=final_alphafold,
    variables=BSC_JOB_VARIABLES + [output, removeExistingResults],
    inputs=[fastaFile],
    outputs=[outputModelsVariable],
)
