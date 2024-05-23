import os
import shutil

from HorusAPI import PluginBlock, PluginVariable, VariableGroup, VariableTypes

# Input variables
pdbFolderVariable = PluginVariable(
    id="pdb_folder",
    name="PDB Folder",
    description="Folder containing the PDB files",
    type=VariableTypes.FOLDER,
)

singlePDBVariable = PluginVariable(
    id="single_pdb",
    name="Single PDB",
    description="Single PDB file to convert",
    type=VariableTypes.FILE,
    allowedValues=["pdb"],
)

structureVariable = PluginVariable(
    id="structure",
    name="Structure",
    description="Structure to convert",
    type=VariableTypes.STRUCTURE,
)

# Internal variables
changeLigandNameVariable = PluginVariable(
    id="change_ligand_name",
    name="Change ligand name",
    description="Change the ligand name inside the PDB. This will replace the chain, residue and atom names with the ligand name (L)",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)

# Output variables
outputVariable = PluginVariable(
    id="output",
    name="MAE folder",
    description="Output folder containing the MAE files",
    type=VariableTypes.FOLDER,
)


def convertPDBToMAE(block: PluginBlock):
    # Test if we have valid glide installation
    command = "echo $SCHRODINGER"
    output = block.remote.remoteCommand(command)

    if output is None or output == "":
        raise Exception(f"No valid Schrodinger installation found on remote {block.remote.name}")
    else:
        print(f"Schrodinger installation found on remote {block.remote.name}: {output}")

    run_command = str(output) + "/run"

    import prepare_proteins

    pdb_folder = block.inputs.get("pdb_folder", None)

    if pdb_folder is None:
        raise Exception("No PDB folder selected")

    if not os.path.isdir(pdb_folder):
        raise Exception(f"Invalid PDB folder: {pdb_folder}")

    # If there are subfolders inside the PDB folder, we need to move the files
    # to the main folder
    hasSubfolders = False
    for subfolder in os.listdir(pdb_folder):
        if os.path.isdir(os.path.join(pdb_folder, subfolder)):
            hasSubfolders = True

    if hasSubfolders:
        # Recursively get all the PDBs inside the folder
        pdb_files = []
        for root, dirs, files in os.walk(pdb_folder):
            for file in files:
                if file.endswith(".pdb"):
                    pdb_files.append(os.path.join(root, file))

        if len(pdb_files) == 0:
            raise Exception(f"No PDB files found in {pdb_folder}")

        # Create a new folder for the PDB files
        folder_name = os.path.basename(pdb_folder)
        pdb_folder = os.path.join(os.getcwd(), f"gathered_pdb_{folder_name}")

        os.makedirs(pdb_folder, exist_ok=True)

        # Copy the PDB files to the new folder
        for pdb_file in pdb_files:
            shutil.copy(pdb_file, pdb_folder)

    models = prepare_proteins.proteinModels(pdb_folder)

    change_ligand_name = block.variables.get("change_ligand_name", False)

    # The first time we run this, the script will be generated but not executed
    # as we don't have the Schrodinger license locally. We need to run it again
    # on the remote once the convert script is generated.
    print("Generating conversion script")
    models.convertLigandPDBtoMae(pdb_folder, change_ligand_name=change_ligand_name)

    mae_folder = os.path.join(os.getcwd(), f"{pdb_folder}_mae")

    # Move the MAE files to the output folder
    os.makedirs(mae_folder, exist_ok=True)

    if block.remote.name != "Local":
        # Upload the folder to the remote
        remote_dir = os.path.join(block.remote.workDir, block.flow.savedID)

        # Create the remote folder, delete it first if it exists
        block.remote.remoteCommand(f"rm -rf {remote_dir}")
        block.remote.remoteCommand(f"mkdir -p {remote_dir}")

        print(f"Uploading {pdb_folder} to {remote_dir}")
        final_remote_dir = block.remote.sendData(pdb_folder, remote_dir)

        # Mock the os.system call
        old_system = os.system

        def mockSystem(command):
            command = f"cd {final_remote_dir} && {command.replace('run', run_command)}"

            block.remote.remoteCommand(command)

        os.system = mockSystem
        print("Running conversion script on remote")
        models.convertLigandPDBtoMae(pdb_folder, change_ligand_name=change_ligand_name)

        # Restore the os.system call
        os.system = old_system

        # Download the output
        print("Downloading MAE files")
        downloaded_path = block.remote.getData(final_remote_dir, os.getcwd())

        # Remove the remote folder
        block.remote.remoteCommand(f"rm -rf {remote_dir}")

        for model in os.listdir(downloaded_path):
            if model.endswith(".mae"):
                os.rename(os.path.join(pdb_folder, model), os.path.join(mae_folder, model))

    elif block.remote.name == "Local":
        for model in os.listdir(pdb_folder):
            if model.endswith(".mae"):
                # Move the MAE files to the output folder
                shutil.move(os.path.join(pdb_folder, model), os.path.join(mae_folder, model))

    print(
        f"Sucessfully converted PDB files to MAE. Files converted: {len(os.listdir(mae_folder))}"
    )

    block.setOutput("output", mae_folder)


convertPDBToMAEBlock = PluginBlock(
    name="PDB to MAE",
    description="Convert PDB files to MAE for Glide",
    inputGroups=[
        VariableGroup(
            id=structureVariable.id,
            name=structureVariable.name,
            description=structureVariable.description,
            variables=[structureVariable],
        ),
        VariableGroup(
            id=singlePDBVariable.id,
            name=singlePDBVariable.name,
            description=singlePDBVariable.description,
            variables=[singlePDBVariable],
        ),
        VariableGroup(
            id=pdbFolderVariable.id,
            name=pdbFolderVariable.name,
            description=pdbFolderVariable.description,
            variables=[pdbFolderVariable],
        ),
    ],
    variables=[changeLigandNameVariable],
    outputs=[outputVariable],
    action=convertPDBToMAE,
)
