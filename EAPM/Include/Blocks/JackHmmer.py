"""
Module containing the JackHmmer block for the EAPM plugin as a nord3 implementation
"""

import os

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
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/shared/public/AlphaFold/uniref90/uniref90.fa",
)


def runJackHmmer(block: SlurmBlock):

    inputfasta = block.inputs.get("input_fasta", None)

    if "nord3" not in block.remote.host:
        raise Exception("This block only works on Nord3.")

    if inputfasta is None:
        raise Exception("No input fasta provided")
    if not os.path.exists(inputfasta):
        raise Exception(f"The input fasta file does not exist: {inputfasta}")

    folderName = block.variables.get("folder_name", "jackHmmer")
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

    output = block.outputs.get("output", "output.hmm")
    sequenceDB = block.variables.get(
        "sequence_db", "/gpfs/projects/shared/public/AlphaFold/uniref90/uniref90.fa"
    )

    jobs = [f"jackhmmer -o {folderName}/{output} {folderName}/{inputfasta} {sequenceDB}"]

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
    from utils import downloadResultsAction

    downloaded_path = downloadResultsAction(block)

    resultsFolder = block.extraData["folder_name"]

    output_hmm = os.path.join(downloaded_path, resultsFolder, "output.hmm")

    block.setOutput(outputVariable.id, output_hmm)


from utils import BSC_JOB_VARIABLES

jackHmmerBlock = SlurmBlock(
    name="HmmAlign",
    initialAction=runJackHmmer,
    finalAction=finalAction,
    description="Align sequences to a profile HMM",
    inputs=[fastaInput],
    variables=BSC_JOB_VARIABLES + [sequenceDBVar, removeExistingResults],
    outputs=[outputVariable],
)
