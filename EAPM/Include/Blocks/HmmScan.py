"""
Module containing the HmmScan block for the EAPM plugin as a nord3 implementation
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

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


def runHmmScan(block: SlurmBlock):
    import os

    input = block.inputs.get("input_fasta", None)

    if "nord3" not in block.remote.host:
        raise Exception("This block only works on Nord3.")

    if input is None:
        raise Exception("No input fasta provided")
    if not os.path.exists(input):
        raise Exception(f"The input fasta file does not exist: {input}")

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
    os.system(f"cp {input} {folderName}")

    hmmDB = block.variables.get("hmm_db", None)
    output = block.outputs.get("output", "output.hmm")

    if block.remote.isLocal:
        hmmerExecutable = block.config.get("hmmer_path", "hmmer") + "/hmmscan"
    else:
        hmmerExecutable = "hmmscan"

    jobs = [f"{hmmerExecutable} {hmmDB} {folderName}/{input} -o {folderName}/{output}"]

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

    output_search = os.path.join(downloaded_path, resultsFolder, "output_models")

    block.setOutput(outputVariable.id, output_search)


from utils import BSC_JOB_VARIABLES

hmmScanBlock = SlurmBlock(
    name="HmmScan",
    id="HmmScan",
    initialAction=runHmmScan,
    finalAction=finalAction,
    description="Search sequence(s) against a profile database",
    inputs=[fastaInput],
    variables=BSC_JOB_VARIABLES + [hmmDBVar, removeExistingResults],
    outputs=[outputVariable],
)
