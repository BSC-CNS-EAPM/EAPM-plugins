"""
Module containing the HmmScan block for the EAPM plugin as a nord3 implementation
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes
from utils import BSC_JOB_VARIABLES

# ==========================#
# Variable inputs
# ==========================#
fastaInput = PluginVariable(
    id="input_fasta",
    name="Fasta input",
    description="The input fast",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)


# ==========================#
# Variable outputs
# ==========================#
outputVariable = PluginVariable(
    id="output",
    name="Output File",
    description="Output of the HmmScan block",
    type=VariableTypes.FILE,
    defaultValue="output.hmm",
)

# ==========================#
# Variable
# ==========================#
removeExistingResults = PluginVariable(
    name="Remove existing results",
    id="remove_existing_results",
    description="Remove existing results",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
hmmDBVar = PluginVariable(
    id="hmm_db",
    name="Hmm DB",
    description="The sequence database to search",
    type=VariableTypes.STRING,
    defaultValue=" ",
)
folderName = PluginVariable(
    name="Folder name",
    id="folder_name",
    description="Folder name for the results",
    type=VariableTypes.STRING,
    defaultValue="hmmScan",
)


def run_hmmscan(block: SlurmBlock):
    """
    Runs the HMMScan calculation for a given SlurmBlock.

    Args:
        block (SlurmBlock): The SlurmBlock object containing the necessary inputs and configuration.

    Raises:
        Exception: If the block is not running on Nord3.
        Exception: If no input fasta file is provided.
        Exception: If the input fasta file does not exist.
        Exception: If the specified folder already exists and remove_existing_results is False.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

    input_fasta = block.inputs.get(fastaInput.id, None)

    if "nord3" not in block.remote.host:
        raise ValueError("This block only works on Nord3.")

    if input is None:
        raise ValueError("No input fasta provided")
    if not os.path.exists(input):
        raise ValueError(f"The input fasta file does not exist: {input_fasta}")

    folder_name = block.variables.get(folderName.id, "hmmScan")
    block.extraData["folder_name"] = folder_name
    remove_existing = block.variables.get(removeExistingResults.id, False)

    # If folder already exists, raise exception
    if remove_existing and os.path.exists(folder_name):
        os.system("rm -rf " + folder_name)

    if not remove_existing and os.path.exists(folder_name):
        raise ValueError(
            f"The folder {folder_name} already exists. Please, choose another name or remove it."
        )

    # Create an copy the inputs
    os.makedirs(folder_name, exist_ok=True)
    os.system(f"cp {input_fasta} {folder_name}")

    hmm_db = block.variables.get(hmmDBVar.id, None)
    output = block.outputs.get(outputVariable.id, "output.hmm")

    if block.remote.isLocal:
        hmmer_executable = block.config.get("hmmer_path", "hmmer") + "/hmmscan"
    else:
        hmmer_executable = "hmmscan"

    jobs = [f"{hmmer_executable} {hmm_db} {folder_name}/{input_fasta} -o {folder_name}/{output}"]

    launchCalculationAction(
        block,
        jobs,
        program="hmmer",
        uploadFolders=[
            folder_name,
        ],
    )


def final_action(block: SlurmBlock):
    """
    Perform the final action for the given SlurmBlock.

    Args:
        block (SlurmBlock): The SlurmBlock object.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import downloadResultsAction

    # pylint: enable=import-outside-toplevel

    downloaded_path = downloadResultsAction(block)

    results_folder = block.extraData["folder_name"]

    output_search = os.path.join(downloaded_path, results_folder, "output_models")

    block.setOutput(outputVariable.id, output_search)


hmmScanBlock = SlurmBlock(
    name="HmmScan",
    id="HmmScan",
    initialAction=run_hmmscan,
    finalAction=final_action,
    description="Search sequence(s) against a profile database",
    inputs=[fastaInput],
    variables=BSC_JOB_VARIABLES + [hmmDBVar, removeExistingResults, folderName],
    outputs=[outputVariable],
)
