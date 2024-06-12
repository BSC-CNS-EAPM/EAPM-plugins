"""
Module containing the PDB2MAE block for the EAPM plugin
"""

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
    description="Change the ligand name inside the PDB. "
    "This will replace the chain, residue and atom names with the ligand name (L)",
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


def convert_pdb_2_mae(block: PluginBlock):
    """
    Converts PDB files to MAE format.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.

    Raises:
        Exception: If no valid Schrödinger installation is found on the remote.
        Exception: If no PDB file is selected.
        Exception: If the selected PDB file is invalid.
        Exception: If no PDB folder is selected.
        Exception: If the selected PDB folder is invalid.
        Exception: If no PDB files are found in the selected folder.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os
    import shutil

    import prepare_proteins

    # pylint: enable=import-outside-toplevel
    # Test if we have valid glide installation
    command = "echo $SCHRODINGER"
    output = block.remote.remoteCommand(command)

    if output is None or output == "":
        raise ValueError(f"No valid Schrödinger installation found on remote {block.remote.name}")
    else:
        print(f"Schrödinger installation found on remote {block.remote.name}: {output}")

    run_command = str(output) + "/run"

    if block.selectedInputGroup == singlePDBVariable.id:
        pdb_file = block.inputs.get("single_pdb", None)

        if pdb_file is None:
            raise ValueError("No PDB file selected")

        if not os.path.isfile(pdb_file):
            raise ValueError(f"Invalid PDB file: {pdb_file}")

        if os.path.exists("tmp_ligand"):
            shutil.rmtree("tmp_ligand")
        os.mkdir("tmp_ligand")
        shutil.copy(pdb_file, "tmp_ligand")

        pdb_folder = os.path.join(os.getcwd(), "tmp_ligand")
    else:
        pdb_folder = block.inputs.get("pdb_folder", None)

    if pdb_folder is None:
        raise ValueError("No PDB folder selected")

    if not os.path.isdir(pdb_folder):
        raise ValueError(f"Invalid PDB folder: {pdb_folder}")

    # If there are subfolders inside the PDB folder, we need to move the files
    # to the main folder
    has_subfolders = False
    for subfolder in os.listdir(pdb_folder):
        if os.path.isdir(os.path.join(pdb_folder, subfolder)):
            has_subfolders = True

    if has_subfolders:
        # Recursively get all the PDBs inside the folder
        pdb_files = []
        for root, _, files in os.walk(pdb_folder):
            for file in files:
                if file.endswith(".pdb"):
                    pdb_files.append(os.path.join(root, file))

        if len(pdb_files) == 0:
            raise ValueError(f"No PDB files found in {pdb_folder}")

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
    # as we don't have the Schrödinger license locally. We need to run it again
    # on the remote once the convert script is generated.
    print("Generating conversion script")
    models.convertLigandPDBtoMae(pdb_folder, change_ligand_name=change_ligand_name)

    mae_folder = os.path.join(os.getcwd(), f"{pdb_folder}_mae")

    # Move the MAE files to the output folder
    os.makedirs(mae_folder, exist_ok=True)

    for model in os.listdir(pdb_folder):
        if model.endswith(".mae"):
            # Move the MAE files to the output folder
            shutil.move(os.path.join(pdb_folder, model), os.path.join(mae_folder, model))

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

        def mock_system(command):
            command = f"cd {final_remote_dir} && {command.replace('run', run_command)}"

            block.remote.remoteCommand(command)

        os.system = mock_system
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
        f"Successfully converted PDB files to MAE. Files converted: {len(os.listdir(mae_folder))}"
    )

    block.setOutput("output", mae_folder)


convertPDBToMAEBlock = PluginBlock(
    name="PDB to MAE",
    id="PDBToMAE",
    description="Convert PDB files to MAE for Glide",
    inputGroups=[
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
        VariableGroup(
            id=structureVariable.id,
            name=structureVariable.name,
            description=structureVariable.description,
            variables=[structureVariable],
        ),
    ],
    variables=[changeLigandNameVariable],
    outputs=[outputVariable],
    action=convert_pdb_2_mae,
)
