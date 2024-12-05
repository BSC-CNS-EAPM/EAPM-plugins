"""
Module containing the PrepWizard block for the EAPM plugin
"""

from HorusAPI import PluginVariable, SlurmBlock, VariableGroup, VariableTypes
from utils import BSC_JOB_VARIABLES

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
    id="folder_name_variable",
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
def prepwizard_action(block: SlurmBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    # pylint: disable=import-outside-toplevel
    import os
    import time
    import traceback

    import prepare_proteins
    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

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
        raise ValueError("No input selected")

    # Get prepWizard variables
    folder_name = block.variables.get(folderNameVariable.id, "prepared_proteins")
    if os.path.exists(folder_name):
        folder_name = folder_name + "_" + str(time.time())
    block.extraData[folderNameVariable.id] = folder_name
    print("Folder name: ", folder_name)
    ph = int(block.variables.get(phPW.id, 7))
    epik_ph = block.variables.get(epikPHPW.id, False)
    sample_water = block.variables.get(sampleWaterPW.id, False)
    remove_hydrogens = block.variables.get(removeHydrogensPW.id, False)
    del_water_hbond_cutoff = block.variables.get(delWaterHbondCutOffPW.id, False)
    fill_loops = block.variables.get(fillLoopsPW.id, False)
    protonation_states = block.variables.get(protonationStatesPW.id, None)
    noepik = block.variables.get(noepikPW.id, False)
    no_prot_assign = block.variables.get(noProtAssignPW.id, False)

    print("Loading pdbs files...")

    models = prepare_proteins.proteinModels(input_folder)

    print("Setting up PrepWizard Optimitzations...")

    folder_name_wizard = folder_name + "_wizard"

    # if block.remote.name.lower() == "local":
    #     prime = False
    # else:
    #     prime = True

    try:
        jobs = models.setUpPrepwizardOptimization(
            prepare_folder=folder_name_wizard,
            pH=ph,
            epik_pH=epik_ph,
            samplewater=sample_water,
            remove_hydrogens=remove_hydrogens,
            delwater_hbond_cutoff=del_water_hbond_cutoff,
            fill_loops=fill_loops,
            protonation_states=protonation_states,
            noepik=noepik,
            noprotassign=no_prot_assign,
            # prime=prime,
        )
    except ValueError as exc:

        trace = traceback.format_exc()

        raise ValueError(f"Error setting up PrepWizard optimization: {trace}") from exc

    print("Jobs ready to be run.")

    launchCalculationAction(block, jobs, "schrodinger", [folder_name_wizard])


def final_prepwizard(block: SlurmBlock):
    """
    Downloads the results of the PrepWizard job and moves the prepared
    proteins to the output folder.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the PrepWizard job.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os
    import shutil

    from utils import downloadResultsAction

    # pylint: enable=import-outside-toplevel

    downloadResultsAction(block)

    folder_name = block.extraData[folderNameVariable.id]

    # Create the output folder containing the prepared proteins
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

    # Move the prepared proteins to the output folder
    for model in os.listdir(os.path.join(folder_name + "_wizard", "output_models")):
        for file in os.listdir(
            os.path.join(folder_name + "_wizard", "output_models", model)
        ):
            if file.endswith(".pdb"):
                final_path = os.path.join(folder_name, file)
                pdb_path = os.path.join(
                    folder_name + "_wizard", "output_models", model, file
                )
                shutil.copyfile(pdb_path, final_path)

    block.setOutput(outputPDB.id, final_path)
    block.setOutput(outputPW.id, folder_name)


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
    id="PrepWizard",
    description="Run Preparation Wizard optimization. (For AMD cluster, workstations and local)",
    initialAction=prepwizard_action,
    finalAction=final_prepwizard,
    variables=block_variables,
    inputGroups=[folderVariableGroup, fileVariableGroup],
    outputs=[outputPDB, outputPW],
)
