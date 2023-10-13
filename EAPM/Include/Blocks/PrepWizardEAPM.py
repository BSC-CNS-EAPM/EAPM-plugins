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
def initialPrepWizard(block: SlurmBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
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
        print(f"{os.path.join(os.getcwd(), folderName)} : {os.path.join(simRemoteDir)}")
        block.remote.sendData(os.path.join(os.getcwd(), folderName), os.path.join(simRemoteDir))
        scriptPath = os.path.join(simRemoteDir, scriptName)
        print(f"{os.path.join(os.getcwd(), scriptName)} : {scriptPath}")
        block.remote.sendData(os.path.join(os.getcwd(), scriptName), scriptPath)

        print("Data sent to the remote.")

        print("Running the simulation...")

        # Run the simulation
        jobID = block.remote.submitJob(scriptPath)
        print(f"Simulation running with job ID {jobID}. Waiting for it to finish...")
    # * Local
    else:
        print("Running the simulation locally...")

        # Run the simulation
        result = subprocess.run(["bash", scriptName], check=False)
        if result.returncode != 0:
            raise Exception("Error running the local simulation.")


# Block's final action
def finalPrepWizard(block: SlurmBlock):
    """
    Final action of the block. It downloads the results from the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    simRemoteDir = os.path.join(block.remote.workDir, folder)

    print("Alphafold calculation finished, downloading results...")

    destPath = os.path.join(os.getcwd(), folder)

    # Transfer the results from the remote
    block.remote.getData(simRemoteDir, destPath)

    print(f"Results transferred to the local machine at: {destPath}")

    print("Setting output of block to the results directory...")
    # Set the output
    folderName = block.variables.get("folder_name", "alphafold")
    block.setOutput("path", os.path.join(destPath, folderName))


prepWizardBlock = SlurmBlock(
    name="PrepWizard",
    description="Run Preparation Wizard.",
    initialAction=initialPrepWizard,
    finalAction=finalPrepWizard,
    variables=[partitionPW, folderNamePW, cpusPW, scriptNamePW],
    inputs=[fasta_fileAF],
    outputs=[outputAF],
)
