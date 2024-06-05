"""
Module containing the HmmBuild block for the EAPM plugin as a nord3 implementation
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

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


def runHmmBuild(block: SlurmBlock):

    import os

    input = block.inputs.get("input_msa", None)

    if "nord3" not in block.remote.host:
        raise Exception("This block only works on Nord3.")

    if input is None:
        raise Exception("No input msa provided")
    if not os.path.exists(input):
        raise Exception(f"The input msa file does not exist: {input}")

    folderName = block.variables.get("folder_name", "hmmBuild")
    block.extraData["folder_name"] = folderName
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

    # Create an copy the inputs
    os.makedirs(folderName, exist_ok=True)
    os.system(f"cp {input} {folderName}")

    output = block.outputs.get("output", "output.hmm")

    if block.remote.isLocal:
        hmmerExecutable = block.config.get("hmmer_path", "hmmer") + "/hmmbuild"
    else:
        hmmerExecutable = "hmmbuild"

    jobs = [f"{hmmerExecutable} {folderName}/{output} {folderName}/{input}"]

    from utils import launchCalculationAction

    launchCalculationAction(
        block,
        jobs,
        program="hmmer",
        uploadFolders=[
            folderName,
        ],
    )


def finalAction(block: SlurmBlock):
    import os

    from utils import downloadResultsAction

    downloaded_path = downloadResultsAction(block)

    resultsFolder = block.extraData["folder_name"]

    output_hmm = os.path.join(downloaded_path, resultsFolder, "output.hmm")

    block.setOutput(outputVariable.id, output_hmm)


from utils import BSC_JOB_VARIABLES

hmmBuildBlock = SlurmBlock(
    name="HmmBuild",
    initialAction=runHmmBuild,
    finalAction=finalAction,
    description="Creates a hmm from a msa.",
    inputs=[msaInput],
    variables=BSC_JOB_VARIABLES + [removeExistingResults],
    outputs=[outputVariable],
)
