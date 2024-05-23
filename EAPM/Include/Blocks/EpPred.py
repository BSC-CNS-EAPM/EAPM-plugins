import datetime
import os
import subprocess

from HorusAPI import PluginVariable, SlurmBlock, VariableList, VariableTypes

# TODO Making the block to work in marenostrum, if not, will work in local.
# TODO Add to documentation
# For the mn execution set default paths

# ==========================#
# Variable inputs
# ==========================#
inputFasta = PluginVariable(
    name="Input fasta",
    id="input_fasta",
    description="The input fasta file. (-i)",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# ==========================#
# Variable outputs
# ==========================#
outputEppred = PluginVariable(
    name="EP-Pred output",
    id="path",
    description="The folder containing the results.",
    type=VariableTypes.FOLDER,
)

##############################
#       Other variables      #
##############################
removeExistingResults = PluginVariable(
    name="Remove existing results",
    id="remove_existing_results",
    description="Remove existing results",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
pssmDir = PluginVariable(
    name="PSSM directory",
    id="pssm_dir",
    description="The directory containing the PSSM files. (-p)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
fastadir = PluginVariable(
    name="Fasta directory",
    id="fasta_dir",
    description="The directory containing the fasta files. (-f)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
ifeatureDir = PluginVariable(
    name="Ifeature directory",
    id="ifeature_dir",
    description="The directory containing the ifeature files. (-id)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
possumDir = PluginVariable(
    name="Possum directory",
    id="possum_dir",
    description="The directory containing the possum files. (-Po)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
ifeatureOut = PluginVariable(
    name="Ifeature out",
    id="ifeature_out",
    description="The directory where the ifeature features are. (-io)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
possumOut = PluginVariable(
    name="Possum out",
    id="possum_out",
    description="The directory for the possum extractions. (-po)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
filteredOut = PluginVariable(
    name="Filtered output",
    id="filtered_out",
    description="The directory for the filtered features. (-fo)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
dbinp = PluginVariable(
    name="Database input",
    id="dbinp",
    description="The path to the fasta files to create the database. (-di)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
dbout = PluginVariable(
    name="Database output",
    id="dbout",
    description="The path and name of the created database. (-do)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
numThread = PluginVariable(
    name="Number of threads",
    id="num_thread",
    description="The number of threads to use for the generation of pssm profiles and feature extraction. (-n)",
    type=VariableTypes.INTEGER,
    defaultValue=5,
)
resDir = PluginVariable(
    name="Result directory",
    id="res_dir",
    description="The name for the folder where to store the prediction results. (-rs)",
    type=VariableTypes.FOLDER,
    defaultValue="results",
)
numSimilarSamples = PluginVariable(
    name="Number of similar samples",
    id="num_similar_samples",
    description="The number of similar training samples to filter the predictions. (-nss)",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
restart = PluginVariable(
    name="Restart",
    id="restart",
    description="From which part of the process to restart with. (-re)",
    type=VariableTypes.STRING_LIST,
    allowedValues=["feature", "predict"],
)
filterOnly = PluginVariable(
    name="Filter only",
    id="filter_only",
    description="True if you already have the features. (-on)",
    type=VariableTypes.BOOLEAN,
    defaultValue=None,
)
extractionRestart = PluginVariable(
    name="Extraction restart",
    id="extraction_restart",
    description="The file to restart the extraction with. (-er)",
    type=VariableTypes.FILE,
)
long = PluginVariable(
    name="Long",
    id="long",
    description="True when restarting from the long commands. (-lg)",
    type=VariableTypes.BOOLEAN,
    defaultValue=None,
)
run = PluginVariable(
    name="Run",
    id="run",
    description="Run possum or ifeature extraction (-r)",
    type=VariableTypes.STRING,
    defaultValue=None,
    allowedValues=["possum", "ifeature", "both"],
)
start = PluginVariable(
    name="Start",
    id="start",
    description="The starting number. (-st)",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
end = PluginVariable(
    name="End",
    id="end",
    description="The ending number, not included. (-en)",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
sbatchPath = PluginVariable(
    name="Sbatch path",
    id="sbatch_path",
    description="The folder to keep the run files for generating pssm. (-sp)",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)
value = PluginVariable(
    name="Value",
    id="value",
    description="The voting threshold to be considered positive. (-v)",
    type=VariableTypes.NUMBER_LIST,
    defaultValue=None,
    allowedValues=[1, 0.8, 0.5],
)
iterations = PluginVariable(
    name="Iterations",
    id="iterations",
    description="The number of iterations in the PSIBlast. (-iter)",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)


def runEppred(block: SlurmBlock):

    inputfasta = block.inputs.get("input_fasta", None)

    if inputfasta is None:
        raise Exception("No input fasta provided")
    if not os.path.exists(inputfasta):
        raise Exception(f"The input fasta file does not exist: {inputfasta}")

    command = "python -m ep_pred.Launch "
    command += f"-i {inputfasta} "
    command += f"-n {block.variables.get('num_thread', 5)} "
    command += f"-re {block.variables.get('restart', 'predict')} "

    pssm_dir = block.variables.get("pssm_dir", None)
    if pssm_dir is not None:
        command += f"-p {pssm_dir} "
    fasta_dir = block.variables.get("fasta_dir", None)
    if fasta_dir is not None:
        command += f"-f {fasta_dir} "
    ifeature_dir = block.variables.get("ifeature_dir", None)
    if ifeature_dir is not None:
        command += f"-id {ifeature_dir} "
    possum_dir = block.variables.get("possum_dir", None)
    if possum_dir is not None:
        command += f"-Po {possum_dir} "
    ifeature_out = block.variables.get("ifeature_out", None)
    if ifeature_out is not None:
        command += f"-io {ifeature_out} "
    possum_out = block.variables.get("possum_out", None)
    if possum_out is not None:
        command += f"-po {possum_out} "
    filtered_out = block.variables.get("filtered_out", None)
    if filtered_out is not None:
        command += f"-fo {filtered_out} "
    dbinp = block.variables.get("dbinp", None)
    if dbinp is not None:
        command += f"-di {dbinp} "
    dbout = block.variables.get("dbout", None)
    if dbout is not None:
        command += f"-do {dbout} "
    res_dir = block.variables.get("res_dir", None)
    if res_dir is not None:
        command += f"-rs {res_dir} "
    num_similar_samples = block.variables.get("num_similar_samples", None)
    if num_similar_samples is not None:
        command += f"-nss {num_similar_samples} "
    restart = block.variables.get("restart", "feature")
    if restart is not None:
        command += f"-re {restart} "
    filter_only = block.variables.get("filter_only", None)
    if filter_only is not None:
        command += f"-on {filter_only} "
    extraction_restart = block.variables.get("extraction_restart", None)
    if extraction_restart is not None:
        command += f"-er {extraction_restart} "
    long = block.variables.get("long", None)
    if long is not None:
        command += f"-lg {long} "
    run = block.variables.get("run", None)
    if run is not None:
        command += f"-r {run} "
    start = block.variables.get("start", None)
    if start is not None:
        command += f"-st {start} "
    end = block.variables.get("end", None)
    if end is not None:
        command += f"-en {end} "
    sbatch_path = block.variables.get("sbatch_path", None)
    if sbatch_path is not None:
        command += f"-sp {sbatch_path} "
    value = block.variables.get("value", None)
    if value is not None:
        command += f"-v {value} "
    iterations = block.variables.get("iterations", None)
    if iterations is not None:
        command += f"-iter {iterations} "

    jobs = [command]

    folderName = block.variables.get("folder_name", "epPred")
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

    from utils import launchCalculationAction

    launchCalculationAction(
        block,
        jobs,
        program="epPred",
        uploadFolders=[
            folderName,
        ],
    )


def finalAction(block: SlurmBlock):
    pass


from utils import BSC_JOB_VARIABLES

epPredBlock = SlurmBlock(
    name="Ep-pred",
    initialAction=runEppred,
    finalAction=finalAction,
    description="A machine learning program to predict promiscuity of esterases.",
    inputs=[inputFasta],
    variables=BSC_JOB_VARIABLES
    + [
        removeExistingResults,
        pssmDir,
        fastadir,
        ifeatureDir,
        possumDir,
        ifeatureOut,
        possumOut,
        filteredOut,
        dbinp,
        dbout,
        numThread,
        resDir,
        numSimilarSamples,
        restart,
        filterOnly,
        extractionRestart,
        long,
        run,
        start,
        end,
        sbatchPath,
        value,
        iterations,
    ],
    outputs=[outputEppred],
)
