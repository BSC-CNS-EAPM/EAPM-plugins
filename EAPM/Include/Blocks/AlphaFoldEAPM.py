"""
Module containing the AlphaFold block for the EAPM plugin
"""

import datetime
import os
import shutil
import subprocess
import tarfile
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
# Variable outputs
# ==========================#
outputAF = PluginVariable(
    name="Alphafold folder output",
    id="path",
    description="The folder containing the results.",
    type=VariableTypes.FOLDER,
)
outputPDBAF = PluginVariable(
    name="Alphafold pdb output",
    id="out_pdb",
    description="First PDB of the Alphafold.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["pdb"]
)

##############################
#       Other variables      #
##############################
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
confidenceThresholdAF = PluginVariable(
    name="Confidence threshold",
    id="confidence_threshold",
    description="Threshold confidence indicates the maximum confidence score at which to stop the trimming of terminal regions.",
    type=VariableTypes.FLOAT,
    defaultValue=90.0,
)

folder = f"AF_{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}"


# Alphafold action block
def initialAlphafold(block: SlurmBlock):
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
        block.remote.sendData(os.path.join(os.getcwd(), folderName), os.path.join(simRemoteDir))
        scriptPath = os.path.join(simRemoteDir, scriptName)
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
def finalAlphafold(block: SlurmBlock):
    """
    Final action of the block. It downloads the results from the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    folderName = block.variables.get("folder_name", "alphafold")

    if block.remote.name != "local":
        cluster = block.remote.host
    else:
        cluster = "local"

    if cluster != "local":
        block.setOutput("path", folderName)

        simRemoteDir = os.path.join(block.remote.workDir, folder)  # AF_********-*******/

        print("Alphafold calculation finished, downloading results...")

        destPath = os.path.join(os.getcwd(), folder)

        # Transfer the results from the remote
        block.remote.remoteCommand(
            f"cd {block.remote.workDir} && tar czvf {simRemoteDir}.tar.gz {folder}/{folderName}/output_models/*/ranked_0.pdb"
        )

        block.remote.getData(f"{simRemoteDir}.tar.gz", f"{destPath}.tar.gz")

        print(f"Results transferred to the local machine at: {destPath}.tar.gz")

        # Unzip the local file
        with tarfile.open(f"{destPath}.tar.gz", "r:gz") as tar:
            tar.extractall()

        # Delete the tar file
        os.remove(f"{destPath}.tar.gz")

        # Final tweaks to the AF results for a better output
        # Create a structures folder if it does not exists
        if not os.path.exists(f"{os.getcwd()}/models"):
            os.mkdir(f"{os.getcwd()}/models")

        print("Setting up folder...")
        # Copy each alphafold output model (from a specific rank) into the models folder
        rank = 0
        for model in os.listdir(f"{destPath}/alphafold/output_models/"):
            if os.path.exists(
                f"{destPath}/alphafold/output_models/" + model + "/ranked_" + str(rank) + ".pdb"
            ):
                shutil.copyfile(
                    f"{destPath}/alphafold/output_models/"
                    + model
                    + "/ranked_"
                    + str(rank)
                    + ".pdb",
                    f"{os.getcwd()}/models/" + model + ".pdb",
                )

        print("Loading prepare_proteins")
        import prepare_proteins

        print("Loading models")
        models = prepare_proteins.proteinModels(os.path.join(os.getcwd(), "models"))
        print("Trimming the models...")
        confidenceThreshold = block.variables.get("confidence_threshold", None)
        models.removeTerminiByConfidenceScore(confidenceThreshold)
        print("Saving trimmed models")
        models.saveModels("trimmed_models")

        print("Setting output of block to the results directory...")
        
        outPdb = None
        for file in os.listdir(os.path.join(os.getcwd(), "trimmed_models")):
            if file.endswith(".pdb"):
                outPdb = os.path.join(os.getcwd(), "trimmed_models", file)
                break

        # Set the output
        block.setOutput("path", os.path.join(os.getcwd(), "trimmed_models"))
        block.setOutput("out_pdb", outPdb)
    # * Local
    else:
        # Final tweaks to the AF results for a better output
        # Create a structures folder if it does not exists
        if not os.path.exists(f"{os.getcwd()}/models"):
            os.mkdir(f"{os.getcwd()}/models")

        print("Setting up folder...")
        # Copy each alphafold output model (from a specific rank) into the models folder
        rank = 0
        for model in os.listdir("/alphafold/output_models/"):
            if os.path.exists(
                "/alphafold/output_models/" + model + "/ranked_" + str(rank) + ".pdb"
            ):
                shutil.copyfile(
                    "/alphafold/output_models/" + model + "/ranked_" + str(rank) + ".pdb",
                    f"{os.getcwd()}/models/" + model + ".pdb",
                )

        print("Loading prepare_proteins")
        import prepare_proteins

        print("Loading models")
        models = prepare_proteins.proteinModels(os.path.join(os.getcwd(), "models"))
        print("Trimming the models...")
        confidenceThreshold = block.variables.get("confidence_threshold", None)
        models.removeTerminiByConfidenceScore(confidenceThreshold)
        print("Saving trimmed models")
        models.saveModels("trimmed_models")

        print("Setting output of block to the results directory...")

        # Set the output
        block.setOutput("path", os.path.join(os.getcwd(), "trimmed_models"))


alphafoldBlock = SlurmBlock(
    name="Alphafold",
    description="Run Alphafold. (For cte_power, marenostrum and minotauro clusters or local)",
    initialAction=initialAlphafold,
    finalAction=finalAlphafold,
    variables=[partitionAF, folderNameAF, cpusAF, scriptNameAF, confidenceThresholdAF],
    inputs=[fasta_fileAF],
    outputs=[outputAF, outputPDBAF],
)
