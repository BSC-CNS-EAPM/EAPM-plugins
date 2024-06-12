"""
Module containing the HmmSearch block for the EAPM plugin as a nord3 implementation
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


# ==========================#
# Variable outputs
# ==========================#
outputVariable = PluginVariable(
    id="output",
    name="Output File",
    description="Output of the HmmSearch block",
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
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/shared/public/AlphaFold/uniref90/uniref90.fa",
)
evalueVar = PluginVariable(
    id="hmmsearch_evalue",
    name="HmmSearch evalue",
    description="The evalue to use",
    type=VariableTypes.FLOAT,
    defaultValue=0.001,
)
folderName = PluginVariable(
    name="Folder name",
    id="folder_name",
    description="Folder name for the results",
    type=VariableTypes.STRING,
    defaultValue="hmmScan",
)


def run_hmmsearch(block: SlurmBlock):
    """
    Run the hmmsearch program to search for sequence homologs using a hidden Markov model.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the current block.

    Raises:
        Exception: If the block is not running on Nord3.
        Exception: If no input hmm file is provided.
        Exception: If the input hmm file does not exist.
        Exception: If the specified folder already exists and remove_existing_results is False.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

    input_hmm = block.inputs.get(hmmInput.id, None)

    if "nord3" not in block.remote.host:
        raise ValueError("This block only works on Nord3.")

    if input_hmm is None:
        raise ValueError("No input hmm provided")
    if not os.path.exists(input_hmm):
        raise ValueError(f"The input hmm file does not exist: {input_hmm}")

    folder_name = block.variables.get(folderName.id, "hmmSearch")
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
    os.system(f"cp {input_hmm} {folder_name}")

    cpus = block.variables.get("cpus")
    evalue = block.variables.get("hmmsearch_evalue", 0.001)
    output = block.outputs.get("output", "output.hmm")
    sequence_db = block.variables.get(
        sequenceDBVar.id, "/gpfs/projects/shared/public/AlphaFold/uniref90/uniref90.fa"
    )

    if block.remote.isLocal:
        hmmer_executable = block.config.get("hmmer_path", "hmmer") + "/hmmsearch"
    else:
        hmmer_executable = "hmmsearch"

    jobs = [
        f"{hmmer_executable} --cpu {cpus} -E {evalue} {folder_name}/{input_hmm} "
        f"{sequence_db} -o {folder_name}/{output}"
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


hmmSearchBlock = SlurmBlock(
    name="HmmSearch",
    id="HmmSearch",
    initialAction=run_hmmsearch,
    finalAction=final_action,
    description="Searches a sequence database with a given hmm",
    inputs=[hmmInput],
    variables=BSC_JOB_VARIABLES + [sequenceDBVar, removeExistingResults],
    outputs=[outputVariable],
)
