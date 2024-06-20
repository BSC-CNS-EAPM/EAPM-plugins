"""
Module containing the HmmAlign block for the EAPM plugin as a nord3 implementation
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes
from utils import BSC_JOB_VARIABLES

# ==========================#
# Variable inputs
# ==========================#
hmmInput = PluginVariable(
    id="input_hmm",
    name="Hmm input",
    description="The input hmm",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["hmm"],
)
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
folderName = PluginVariable(
    name="Folder name",
    id="folder_name",
    description="Folder name for the results",
    type=VariableTypes.STRING,
    defaultValue="hmmScan",
)


def run_hmmalign(block: SlurmBlock):
    """
    Run the HmmAlign block.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the block.
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

    input_fasta = block.inputs.get(fastaInput.id, None)
    input_hmm = block.inputs.get(hmmInput.id, None)

    if "nord3" not in block.remote.host:
        raise ValueError("This block only works on Nord3.")

    if input_fasta is None:
        raise ValueError("No input fasta provided")
    if not os.path.exists(input_fasta):
        raise ValueError(f"The input fasta file does not exist: {input_fasta}")

    if input_hmm is None:
        raise ValueError("No input Hmm provided")
    if not os.path.exists(input_hmm):
        raise ValueError(f"The input Hmm file does not exist: {input_hmm}")

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
    os.system(f"cp {input_hmm} {folder_name}")

    if block.remote.isLocal:
        hmmer_executable = block.config.get("hmmer_path", "hmmer") + "/hmmalign"
    else:
        hmmer_executable = "hmmalign"

    jobs = [f"{hmmer_executable} {folder_name}/{input_hmm} {folder_name}/{input_fasta}"]

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
    Perform the final action for the SlurmBlock.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the block.
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import downloadResultsAction

    # pylint: enable=import-outside-toplevel

    downloaded_path = downloadResultsAction(block)

    results_folder = block.extraData["folder_name"]

    output_hmm = os.path.join(downloaded_path, results_folder, "output.hmm")

    block.setOutput(outputVariable.id, output_hmm)


hmmAlignBlock = SlurmBlock(
    name="HmmAlign",
    id="hmmAlign",
    initialAction=run_hmmalign,
    finalAction=final_action,
    description="Align sequences to a profile HMM",
    inputs=[hmmInput, fastaInput],
    variables=BSC_JOB_VARIABLES + [removeExistingResults],
    outputs=[outputVariable],
)
