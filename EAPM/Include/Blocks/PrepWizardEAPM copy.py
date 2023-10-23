"""
Module containing the AlphaFold block for the EAPM plugin
"""

import datetime
import os
import subprocess
from HorusAPI import PluginVariable, PluginBlock, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
inputFolderPW = PluginVariable(
    name="Input Folder",
    id="input_folder",
    description="Folder with the pdbs input files.",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)

# ==========================#
# Variable outputs
# ==========================#
outputPW = PluginVariable(
    name="Prepwizard output",
    id="path",
    description="The folder containing the results.",
    type=VariableTypes.FOLDER,
)

##############################
#       Other variables      #
##############################
partitionPW = PluginVariable(
    name="Partition",
    id="partition",
    description="Partition where to lunch.",
    type=VariableTypes.STRING_LIST,
    defaultValue="bsc_ls",
    allowedValues=["bsc_ls", "debug"],
)
cpusPW = PluginVariable(
    name="CPUs",
    id="cpus",
    description="Number of CPUs to use.",
    type=VariableTypes.INTEGER,
    defaultValue=1,
)
workstationPW = PluginVariable(
    name="Workstation",
    id="workstation",
    description="Workstation where to launch the calculation.",
    type=VariableTypes.STRING_LIST,
    defaultValue=["bubbles", "blossom", "cactus"],
)

##############################
# Block's advanced variables #
##############################
folderNamePW = PluginVariable(
    name="Simulation name",
    id="folder_name",
    description="Name of the simulation folder.",
    type=VariableTypes.STRING,
    defaultValue="prepwizard",
)
scriptNamePW = PluginVariable(
    name="Simulation name",
    id="script_name",
    description="Name of the script.",
    type=VariableTypes.STRING,
    defaultValue="slurm_array.sh",
)
phPW = PluginVariable(
    name="PH",
    id="ph",
    description="PH to use.",
    type=VariableTypes.FLOAT,
    defaultValue=7.0,
)
epikPHPW = PluginVariable(
    name="Epik PH",
    id="epik_ph",
    description="If use epik_ph.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
sampleWaterPW = PluginVariable(
    name="Sample Water",
    id="sample_water",
    description="If sample water.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
removeHydrogensPW = PluginVariable(
    name="Remove Hydrogens",
    id="remove_hydrogens",
    description="If remove hydrogens.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
delWaterHbondCutOffPW = PluginVariable(
    name="Delete Water Hbond Cut Off",
    id="del_water_hbond_cut_off",
    description="Delete water hbond cut off.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
fillLoopsPW = PluginVariable(
    name="Fill Loops",
    id="fill_loops",
    description="If fill loops.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
protonationStatesPW = PluginVariable(
    name="Protonation States",
    id="protonation_states",
    description="If protonation states.",
    type=VariableTypes.LIST,
    defaultValue=None,
)
noepikPW = PluginVariable(
    name="No Epik",
    id="no_epik",
    description="If no epik.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
noProtAssignPW = PluginVariable(
    name="No Prot Assign",
    id="no_prot_assign",
    description="If no prot assign.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)

folder = f"PW_{datetime.datetime.now().strftime('%Y%m%d-%H%M%S')}"


# Alphafold action block
def initialPrepWizard(block: PluginBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    # Dictionary with the workstations IPs
    ips = {"bubbles": "84.88.51.219", "cactus": "84.88.51.217", "blossom": "84.88.51.250"}

    # Loading plugin variables
    inputFolder = block.inputs.get("input_folder", "None")
    if inputFolder == "None":
        raise Exception("No folder provided.")
    partition = block.variables.get("partition", "bsc_ls")
    cpus = block.variables.get("cpus", 1)
    workstationPW = block.variables.get("workstation", "bubbles")

    # Get prepWizard variables
    folderName = block.variables.get("folder_name", "prepwizard")
    scriptName = block.variables.get("script_name", "commands")
    ph = block.variables.get("ph", 7.0)
    epikPH = block.variables.get("epik_ph", False)
    sampleWater = block.variables.get("sample_water", False)
    removeHydrogens = block.variables.get("remove_hydrogens", False)
    delWaterHbondCutOff = block.variables.get("del_water_hbond_cut_off", False)
    fillLoops = block.variables.get("fill_loops", False)
    protonationStates = block.variables.get("protonation_states", None)
    noepik = block.variables.get("no_epik", False)
    noProtAssign = block.variables.get("no_prot_assign", False)

    import prepare_proteins

    print("Loading pdbs files...")

    models = prepare_proteins.proteinModels(inputFolder)

    print("Setting up PrepWizard Optimitzations...")

    jobs = models.setUpPrepwizardOptimization(
        prepare_folder=folder,
        pH=ph,
        epik_pH=epikPH,
        samplewater=sampleWater,
        remove_hydrogens=removeHydrogens,
        delwater_hbond_cutoff=delWaterHbondCutOff,
        fill_loops=fillLoops,
        protonation_states=protonationStates,
        noepik=noepik,
        noprotassign=noProtAssign,
    )

    print("Preparing launch files...")

    import bsc_calculations

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


prepWizardBlock = PluginBlock(
    name="PrepWizard AMD cluster",
    description="Run Preparation Wizard optimization. (For AMD cluster)",
    action=initialPrepWizard,
    variables=[
        partitionPW,
        cpusPW,
        folderNamePW,
        scriptNamePW,
        phPW,
        epikPHPW,
        sampleWaterPW,
        removeHydrogensPW,
        delWaterHbondCutOffPW,
        fillLoopsPW,
        protonationStatesPW,
        noepikPW,
        noProtAssignPW,
    ],
    inputs=[inputFolderPW],
    outputs=[outputPW],
)
