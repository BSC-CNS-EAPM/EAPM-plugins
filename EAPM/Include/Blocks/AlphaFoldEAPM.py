"""
Module containing the AlphaFold block for the NBDSuite plugin
"""

import datetime
import os
import subprocess
from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
fasta_fileAF = PluginVariable(
    name="Fasta file",
    id="fasta_file",
    description="The input fasta file.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# ==========================#
# Other variables
# ==========================#
partitionAF = PluginVariable(
    name="Partition",
    id="partition",
    description="Partition where to lunch.",
    type=VariableTypes.STRING_LIST,
    defaultValue="bsc_ls",
    allowedValues=["bsc_ls", "debug"],
)
# clusterAF = PluginVariable(
#     name="Cluster",
#     id="cluster",
#     description="Cluster where to lunch.",
#     type=VariableTypes.STRING_LIST,
#     defaultValue="minotauro",
#     allowedValues=["cte_power", "marenostrum", "minotauro"],
# )
outputAF = PluginVariable(
    name="Alphafold output",
    id="path",
    description="The folder containing the results.",
    type=VariableTypes.FOLDER,
)
cpusAF = PluginVariable(
    name="CPUs",
    id="cpus",
    description="Number of CPUs to use.",
    type=VariableTypes.INTEGER,
    defaultValue=1,
)

##############################
# Block's advanced variables #
##############################
folderNameAF = PluginVariable(
    name="Simulation name",
    id="folder_name",
    description="Name of the simulation folder.",
    type=VariableTypes.STRING,
    defaultValue="alphafold",
)
scriptNameAF = PluginVariable(
    name="Simulation name",
    id="script_name",
    description="Name of the script.",
    type=VariableTypes.STRING,
    defaultValue="slurm_array.sh",
)

folder = f"AF_{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}"


# Alphafold action block
def initialAlphafold(block: SlurmBlock):
    # Loading plugin variables
    fastaFile = block.inputs.get("fasta_file", "None")
    if fastaFile == "None":
        raise Exception("No fasta file provided.")
    partition = block.variables.get("partition", None)
    folderName = block.variables.get("folder_name", None)
    scriptName = block.variables.get("script_name", None)
    cpus = block.variables.get("cpus", None)

    import prepare_proteins

    print("Loading fasta files...")
    print(fastaFile)
    print("dd")
    sequences = prepare_proteins.sequenceModels(fastaFile)

    print("Setting up AlphaFold...")

    jobs = sequences.setUpAlphaFold(folderName)

    print("Preparing launch files...")

    import bsc_calculations

    if block.remote.name != "local":
        cluster = block.remote.host
    else:
        cluster = "local"
        scriptName = "commands"
    print(f"Cluster: {cluster}")
    ## Define cluster
    # cte_power
    if cluster == "plogin1.bsc.es":
        bsc_calculations.cte_power.jobArrays(
            jobs,
            job_name="AF_sequences",
            partition=partition,
            program="alphafold",
            script_name=scriptName,
            cpus=cpus,
        )
    # marenostrum
    elif cluster == "mn1.bsc.es":
        bsc_calculations.marenostrum.jobArrays(
            jobs,
            job_name="AF_sequences",
            partition=partition,
            program="alphafold",
            script_name=scriptName,
            cpus=cpus,
        )
    # minotauro
    elif cluster == "mt1.bsc.es":
        bsc_calculations.minotauro.jobArrays(
            jobs,
            job_name="AF_sequences",
            partition=partition,
            program="alphafold",
            script_name=scriptName,
            gpus=cpus,
        )
    # local
    elif cluster == "local":
        bsc_calculations.local.parallel(
            jobs,
            cpus=min(40, len(jobs)),
        )
    else:
        raise Exception("Cluster not supported.")

    folderName = block.variables.get("folder_name", "alphafold")

    # Create the simulation folder in the remote

    # * Remote cluster
    if cluster != "local":
        simRemoteDir = os.path.join(block.remote.workDir, folder)
        block.remote.remoteCommand(f"mkdir -p -v {simRemoteDir}")

        print(f"Created simulation folder in the remote at {simRemoteDir}")
        print("Sending data to the remote...")
        # Send the system data to the remote
        block.remote.sendData(folderName, os.path.join(simRemoteDir, folderName))
        scriptPath = os.path.join(simRemoteDir, scriptName)
        block.remote.sendData(scriptName, scriptPath)

        print("Data sent to the remote.")

        print("Running the simulation...")

        # Run the simulation
        # jobID = block.remote.submitJob(scriptPath)
        jobID = 1
        print(f"Simulation running with job ID {jobID}. Waiting for it to finish...")
    # * Local
    else:
        print("Running the simulation locally...")

        # Run the simulation
        result = subprocess.run(["bash", scriptName], check=False)
        if result.returncode != 0:
            raise Exception("Error running the local simulation.")


# Block's final action
def finalAction(block: SlurmBlock):
    folderName = block.variables.get("folder_name", "ifd")
    simRemoteDir = os.path.join(block.remote.workDir, folderName)

    print("Alphafold calculation finished, downloading results...")

    destPath = os.path.join(os.getcwd(), folderName)

    # Transfer the results from the remote
    block.remote.getData(simRemoteDir, destPath)

    print(f"Results transferred to the local machine at: {destPath}")

    print("Setting output of block to the results directory...")
    # Set the output
    block.setOutput("path", destPath)


alphafoldBlock = SlurmBlock(
    name="Alphafold",
    description="Run Alphafold.",
    initialAction=initialAlphafold,
    finalAction=finalAction,
    variables=[partitionAF, folderNameAF, cpusAF, scriptNameAF],
    inputs=[fasta_fileAF],
    outputs=[outputAF],
)
