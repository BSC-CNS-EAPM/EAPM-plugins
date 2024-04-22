"""
Module containing the HmmSearch block for the EAPM plugin as a nord3 implementation
"""

import os

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


# ==========================#
# Variable outputs
# ==========================#
outputVariable = PluginVariable(
    id="output",
    name="Output File",
    description="Output of the HmmSearch block",
    type=VariableTypes.FILE,
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
outputHS = PluginVariable(
    name="HmmSearch simulation folder",
    id="folder_name",
    description="The name of the folder where the simulation will be stored.",
    type=VariableTypes.STRING,
    defaultValue="hmmSearch",
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


def runHmmSearch(block: SlurmBlock):

    input = block.inputs.get("input_hmm", None)

    if "nord3" not in block.remote.host:
        raise Exception("This block only works on Nord3.")

    if input is None:
        raise Exception("No input hmm provided")
    if not os.path.exists(input):
        raise Exception(f"The input hmm file does not exist: {input}")

    folderName = block.variables.get("folder_name", "hmmSearch")
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

    cpus = block.variables.get("cpus")
    evalue = block.variables.get("hmmsearch_evalue", 0.001)
    sequenceDBVar = block.variables.get(
        "sequence_db", "/gpfs/projects/shared/public/AlphaFold/uniref90/uniref90.fa"
    )

    jobs = [f"hmmsearch --cpu {cpus} -E {evalue} {input} {sequenceDBVar}"]

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

    output_search = os.path.join(downloaded_path, resultsFolder, "output_models")

    block.setOutput(outputVariable.id, output_search)


from utils import BSC_JOB_VARIABLES

hmmSearchBlock = SlurmBlock(
    name="HmmSearch",
    initialAction=runHmmSearch,
    finalAction=finalAction,
    description="Searches a sequence database with a given hmm",
    inputs=[hmmInput],
    variables=BSC_JOB_VARIABLES + [outputHS, removeExistingResults],
    outputs=[outputVariable],
)
