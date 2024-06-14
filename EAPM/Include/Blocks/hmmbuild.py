"""
Module containing the HmmBuild block for the EAPM plugin as a nord3 implementation
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes
from utils import BSC_JOB_VARIABLES

# ==========================#
# Variable inputs
# ==========================#
msaInput = PluginVariable(
    id="input_msa",
    name="Input MSA",
    description="The input msa",
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
    description="Output of the HmmBuild block",
    type=VariableTypes.FILE,
    defaultValue="output.hmm",
    allowedValues=["hmm"],
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


def run_hmm_build(block: SlurmBlock):
    """
    Run the HmmBuild block.

    Args:
        block (SlurmBlock): The SlurmBlock object.
    """

    # pylint: disable=import-outside-toplevel
    import os

    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

    input_msa = block.inputs.get(msaInput.id, None)

    if "nord3" not in block.remote.host:
        raise ValueError("This block only works on Nord3.")

    if input_msa is None:
        raise ValueError("No input msa provided")
    if not os.path.exists(input_msa):
        raise ValueError(f"The input msa file does not exist: {input_msa}")

    folder_name = block.variables.get(folderName.id, "hmmBuild")
    block.extraData["folder_name"] = folderName
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
    os.system(f"cp {input} {folder_name}")

    output = block.outputs.get("output", "output.hmm")

    if block.remote.isLocal:
        hmmer_executable = block.config.get("hmmer_path", "hmmer") + "/hmmbuild"
    else:
        hmmer_executable = "hmmbuild"

    jobs = [f"{hmmer_executable} {folder_name}/{output} {folder_name}/{input}"]

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


hmmBuildBlock = SlurmBlock(
    name="HmmBuild",
    id="HmmBuild",
    initialAction=run_hmm_build,
    finalAction=final_action,
    description="Creates a hmm from a msa.",
    inputs=[msaInput],
    variables=BSC_JOB_VARIABLES + [removeExistingResults, folderName],
    outputs=[outputVariable],
)
