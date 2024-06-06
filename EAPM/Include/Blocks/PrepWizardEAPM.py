"""
Module containing the PrepWizard block for the EAPM plugin
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableGroup, VariableTypes

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
inputFilePW = PluginVariable(
    name="Input File",
    id="input_file",
    description="File of the pdb to prepare.",
    type=VariableTypes.FILE,
    allowedValues=["pdb"],
)
folderVariableGroup = VariableGroup(
    id="folder_variable_group",
    name="Folder variable group",
    description="Input folder with the models.",
    variables=[inputFolderPW],
)
fileVariableGroup = VariableGroup(
    id="file_output_variable_group",
    name="PDB file group",
    description="Input PDB file.",
    variables=[inputFilePW],
)


# ==========================#
# Variable outputs
# ==========================#
outputPW = PluginVariable(
    name="Prepared proteins",
    id="prepared_proteins",
    description="Folder containing the prepared proteins.",
    type=VariableTypes.FOLDER,
)
outputPDB = PluginVariable(
    name="Output PDB",
    id="out_pdb",
    description="Last PDB of the Prepwizard.",
    type=VariableTypes.FILE,
    allowedValues=["pdb"],
)

##############################
# Block's advanced variables #
##############################
folderNameVariable = PluginVariable(
    name="Folder Name",
    id="folder_name",
    description="Folder name for the prepared proteins.",
    type=VariableTypes.STRING,
    defaultValue="prepared_proteins",
)


# Variables
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


# PrepWizard action block
def prepWizardAction(block: SlurmBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    import os
    import time

    if block.selectedInputGroup == fileVariableGroup.id:
        input_file = block.inputs.get(inputFilePW.id, None)
        input_folder = "models"
        if os.path.exists(input_folder):
            input_folder = input_folder + "_" + str(time.time())
        os.makedirs(input_folder, exist_ok=True)
        os.system(f"cp {input_file} {input_folder}")
    elif block.selectedInputGroup == folderVariableGroup.id:
        input_folder = block.inputs.get(inputFolderPW.id, None)
    else:
        raise Exception("No input selected")

    # Get prepWizard variables
    folderName = block.variables.get(folderNameVariable.id, "prepared_proteins")
    if os.path.exists(folderName):
        folderName = folderName + "_" + str(time.time())
        block.extraData[folderNameVariable.id] = folderName
    ph = int(block.variables.get(phPW.id, 7))
    epikPH = block.variables.get(epikPHPW.id, False)
    sampleWater = block.variables.get(sampleWaterPW.id, False)
    removeHydrogens = block.variables.get(removeHydrogensPW.id, False)
    delWaterHbondCutOff = block.variables.get(delWaterHbondCutOffPW.id, False)
    fillLoops = block.variables.get(fillLoopsPW.id, False)
    protonationStates = block.variables.get(protonationStatesPW.id, None)
    noepik = block.variables.get(noepikPW.id, False)
    noProtAssign = block.variables.get(noProtAssignPW.id, False)

    import prepare_proteins

    print("Loading pdbs files...")

    models = prepare_proteins.proteinModels(input_folder)

    print("Setting up PrepWizard Optimitzations...")

    folderNameWizard = folderName + "_wizard"

    if block.remote.name.lower() == "local":
        prime = False
    else:
        prime = True

    try:
        jobs = models.setUpPrepwizardOptimization(
            prepare_folder=folderNameWizard,
            pH=ph,
            epik_pH=epikPH,
            samplewater=sampleWater,
            remove_hydrogens=removeHydrogens,
            delwater_hbond_cutoff=delWaterHbondCutOff,
            fill_loops=fillLoops,
            protonation_states=protonationStates,
            noepik=noepik,
            noprotassign=noProtAssign,
            # prime=prime,
        )
    except Exception as exc:
        import traceback

        trace = traceback.format_exc()

        raise Exception(f"Error setting up PrepWizard optimization: {trace}") from exc

    print("Jobs ready to be run.")

    from utils import launchCalculationAction

    launchCalculationAction(block, jobs, "schrodinger", [folderNameWizard])


def downloadPrepWizardResults(block: SlurmBlock):
    import os

    from utils import downloadResultsAction

    downloadResultsAction(block)

    folderName = block.extraData[folderNameVariable.id]

    # Create the output folder containing the prepared proteins
    if not os.path.exists(folderName):
        os.mkdir(folderName)

    import shutil

    # Move the prepared proteins to the output folder
    for model in os.listdir(os.path.join(folderName + "_wizard", "output_models")):
        for file in os.listdir(os.path.join(folderName + "_wizard", "output_models", model)):
            if file.endswith(".pdb"):
                finalPath = os.path.join(folderName, file)
                pdbPath = os.path.join(folderName + "_wizard", "output_models", model, file)
                shutil.copyfile(pdbPath, finalPath)

    block.setOutput(outputPDB.id, finalPath)
    block.setOutput(outputPW.id, folderName)


from utils import BSC_JOB_VARIABLES

block_variables = BSC_JOB_VARIABLES + [
    folderNameVariable,
    phPW,
    epikPHPW,
    sampleWaterPW,
    removeHydrogensPW,
    delWaterHbondCutOffPW,
    fillLoopsPW,
    protonationStatesPW,
    noepikPW,
    noProtAssignPW,
]


prepWizardBlock = SlurmBlock(
    name="PrepWizard",
    description="Run Preparation Wizard optimization. (For AMD cluster, workstations and local)",
    initialAction=prepWizardAction,
    finalAction=downloadPrepWizardResults,
    variables=block_variables,
    inputGroups=[folderVariableGroup, fileVariableGroup],
    outputs=[outputPDB, outputPW],
)
