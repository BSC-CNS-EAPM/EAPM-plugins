"""
Module containing the PrepWizard block for the EAPM plugin
"""

import os
from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

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
    name="Prepared proteins",
    id="prepared_proteins",
    description="Folder containing the prepared proteins.",
    type=VariableTypes.FOLDER,
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
    # Loading plugin variables
    inputFolder = block.inputs.get("input_folder", None)
    if inputFolder is None:
        raise Exception("No input folder provided.")

    # Get prepWizard variables
    folderName = block.variables.get("folder_name", "prepared_proteins")
    ph = int(block.variables.get("ph", 7))
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
            prime=prime,
        )
    except Exception as exc:
        import traceback

        trace = traceback.format_exc()

        raise Exception(f"Error setting up PrepWizard optimization: {trace}") from exc

    print("Jobs ready to be run.")

    from utils import launchCalculationAction

    launchCalculationAction(block, jobs, "schrodinger", [folderNameWizard])


def downloadPrepWizardResults(block: SlurmBlock):
    from utils import downloadResultsAction

    downloadResultsAction(block)

    folderName = block.variables.get("folder_name")

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

    block.setOutput("prepared_proteins", folderName)


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
    inputs=[inputFolderPW],
    outputs=[outputPW],
)
