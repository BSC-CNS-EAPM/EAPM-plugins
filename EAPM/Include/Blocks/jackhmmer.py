"""
Module containing the JackHmmer block for the EAPM plugin as a nord3 implementation
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
    description="Output of the JackHmmer block",
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
sequenceDBVar = PluginVariable(
    id="sequence_db",
    name="Sequence DB",
    description="The sequence database to search",
    type=VariableTypes.FILE,
    defaultValue="/apps/ACC/ALPHAFOLD/SRC/database/Alphafold/uniref90/uniref90.fasta",
)
folderNameVar = PluginVariable(
    id="folder_name",
    name="Folder name",
    description="The folder name",
    type=VariableTypes.STRING,
    defaultValue="jackHmmer",
)


def run_jackhmmer(block: SlurmBlock):
    """
    Runs the jackhmmer program to perform sequence alignment and profile construction.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the current block.

    Raises:
        ValueError: If no input fasta file is provided.
        FileNotFoundError: If the input fasta file does not exist.
        Exception: If the folder already exists and remove_existing_results is set to False.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

    input_fasta = block.inputs.get("input_fasta", None)

    # if "nord3" not in block.remote.host or "glogin" not in block.remote.host:
    #     raise Exception("This block only works on Nord3 or mn.")

    if input_fasta is None:
        raise ValueError("No input fasta provided")
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"The input fasta file does not exist: {input_fasta}")

    folder_name = block.variables.get(folderNameVar.id, "jackHmmer")
    block.extraData["folder_name"] = folder_name
    remove_existing = block.variables.get(removeExistingResults.id, False)

    # If folder already exists, raise exception
    if remove_existing and os.path.exists(folder_name):
        os.system("rm -rf " + folder_name)

    if not remove_existing and os.path.exists(folder_name):
        raise FileExistsError(
            f"The folder {folder_name} already exists. Please, "
            "choose another name or remove it with the RemoveExistingFolder option."
        )

    # Create an copy the inputs
    os.makedirs(folder_name, exist_ok=True)
    os.system(f"cp {input_fasta} {folder_name}")

    input_fasta = os.path.join(folder_name, os.path.basename(input_fasta))

    output = block.outputs.get(outputVariable.id, "output.hmm")
    sequence_db = block.variables.get(
        "sequence_db", "/apps/ACC/ALPHAFOLD/SRC/database/Alphafold/uniref90/uniref90.fasta"
    )
    cpus = block.variables.get("cpus", 1)

    if block.remote.isLocal:
        hmmer_executable = block.config.get("hmmer_path", "hmmer") + "/jackhmmer"
    else:
        hmmer_executable = "jackhmmer"
    jobs = [
        f"{hmmer_executable} -o {folder_name}/{output} --cpu {cpus} {input_fasta} {sequence_db}"
    ]

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
        block (SlurmBlock): The SlurmBlock object representing the block.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import downloadResultsAction

    # pylint: enable=import-outside-toplevel

    downloaded_path = downloadResultsAction(block)

    results_folder = block.extraData["folder_name"]

    output_hmm = os.path.join(downloaded_path, results_folder, "output.hmm")

    block.setOutput(outputVariable.id, output_hmm)


jackHmmerBlock = SlurmBlock(
    name="JackHmmer",
    id="JackHmmer",
    initialAction=run_jackhmmer,
    finalAction=final_action,
    description="Iteratively search a protein sequence against a protein database",
    inputs=[fastaInput],
    variables=BSC_JOB_VARIABLES + [sequenceDBVar, removeExistingResults, folderNameVar],
    outputs=[outputVariable],
)
