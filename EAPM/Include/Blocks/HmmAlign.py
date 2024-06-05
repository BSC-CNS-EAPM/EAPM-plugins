"""
Module containing the HmmAlign block for the EAPM plugin as a nord3 implementation
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

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


def runHmmAlign(block: SlurmBlock):

    import os

    inputfasta = block.inputs.get("input_fasta", None)
    inputhmm = block.inputs.get("input_hmm", None)

    if "nord3" not in block.remote.host:
        raise Exception("This block only works on Nord3.")

    if inputfasta is None:
        raise Exception("No input fasta provided")
    if not os.path.exists(inputfasta):
        raise Exception(f"The input fasta file does not exist: {inputfasta}")

    if inputhmm is None:
        raise Exception("No input Hmm provided")
    if not os.path.exists(inputhmm):
        raise Exception(f"The input Hmm file does not exist: {inputhmm}")

    folderName = block.variables.get("folder_name", "hmmScan")
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
    os.system(f"cp {inputfasta} {folderName}")
    os.system(f"cp {inputhmm} {folderName}")

    if block.remote.isLocal:
        hmmerExecutable = block.config.get("hmmer_path", "hmmer") + "/hmmalign"
    else:
        hmmerExecutable = "hmmalign"

    jobs = [f"{hmmerExecutable} {folderName}/{inputhmm} {folderName}/{inputfasta}"]

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

hmmAlignBlock = SlurmBlock(
    name="HmmAlign",
    initialAction=runHmmAlign,
    finalAction=finalAction,
    description="Align sequences to a profile HMM",
    inputs=[hmmInput, fastaInput],
    variables=BSC_JOB_VARIABLES + [removeExistingResults],
    outputs=[outputVariable],
)
