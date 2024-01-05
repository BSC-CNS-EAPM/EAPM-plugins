import os

from HorusAPI import PluginBlock, VariableTypes, PluginVariable, VariableGroup

# Input variables
modelFolderVariable = PluginVariable(
    id="model_folder",
    name="Model folder",
    description="Folder containing the models",
    type=VariableTypes.FOLDER,
    defaultValue="models",
)

ligandFolderVariable = PluginVariable(
    id="ligand_folder",
    name="Ligand folder",
    description="Folder containing the ligands",
    type=VariableTypes.FOLDER,
    defaultValue="ligands",
)

dockingCenterVariable = PluginVariable(
    id="outer_box",
    name="Outer box",
    description="Select the center of the box for docking. This will define the outer box size too.",
    type=VariableTypes.SPHERE,
)

innerBoxVariable = PluginVariable(
    id="inner_box",
    name="Inner Box",
    description="Select the inner box size for docking. From this variable only the radius will be used.",
    type=VariableTypes.SPHERE,
)

# Output variables
outputJobsVariable = PluginVariable(
    id="output_jobs",
    name="BSC jobs",
    description="Output jobs from the docking calculation",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_jobs"],
)


# Action
def glideDocking(block: PluginBlock):
    import prepare_proteins

    models_folder = block.inputs.get("model_folder")

    if models_folder is None or not os.path.isdir(models_folder):
        raise Exception("No valid models folder selected")

    ligand_folder = block.inputs.get("ligand_folder")

    if ligandFolderVariable is None or not os.path.isdir(ligand_folder):
        raise Exception("No valid ligands folder selected")

    docking_center = block.inputs["outer_box"]
    radius = docking_center["radius"]
    center = docking_center["center"]
    x = center["x"]
    y = center["y"]
    z = center["z"]

    for ligand in os.listdir(ligand_folder):
        if ligand.endswith(".pdb"):
            raise Exception(
                "PDB files are not supported as ligands. Please convert to MAE using the 'PDB to MAE' block"
            )

    for model in os.listdir(models_folder):
        if model.endswith(".pdb"):
            raise Exception(
                "PDB files are not supported as models. Please convert to MAE using the 'PDB to MAE' block"
            )

    # Get the PDB models
    pdb_models_folder = models_folder.replace("_mae", "")

    if not os.path.isdir(pdb_models_folder):
        raise Exception(
            f"PDB models folder ({pdb_models_folder}) not found. Please convert the models to PDB using the 'MAE to PDB' block and keep the original PDB models folder"
        )

    models = prepare_proteins.proteinModels(pdb_models_folder)

    # Create a fake center_atoms variable for the library
    center_atoms = {}
    for model in os.listdir(models_folder):
        if model.endswith(".mae"):
            model_name = model.split(".")[0]
            center_atoms[model_name] = [x, y, z]

    # If the center_atoms is empty, raise an exception
    if len(center_atoms) == 0:
        raise Exception("No valid models found in the models folder")

    # Create a box from the radius
    radius = radius
    outerbox = (radius, radius, radius)

    inner_box_size = block.inputs["inner_box"]["radius"]
    innerbox = (inner_box_size, inner_box_size, inner_box_size)

    jobs = models.setUpDockingGrid(
        "grid",
        center_atoms,
        innerbox=innerbox,
        outerbox=outerbox,
    )  # Set grid calcualtion

    # Copy the models .mae to the grid/input_models folder
    for model in os.listdir(models_folder):
        if model.endswith(".mae"):
            model_path = os.path.join(models_folder, model)
            new_model_path = os.path.join("grid/input_models", model)
            os.system(f"cp {model_path} {new_model_path}")

    # Remove any .pdb in the grid/input_models folder
    os.system("rm grid/input_models/*.pdb")

    if len(jobs) == 0:
        raise Exception("No jobs created. Is Glide correctly installed?")

    print(f"Sucessfully prepared docking grid input. Jobs created: {len(jobs)}")

    output_jobs = {
        "program": "glide",
        "jobs": jobs,
        "results_data": {
            "ligand_folder": ligand_folder,
            "model_folder": models_folder,
        },
    }

    block.setOutput("output_jobs", output_jobs)


setupDockingGrid = PluginBlock(
    name="Setup Docking Grid",
    description="Prepare the files for a GLIDE docking calculation",
    action=glideDocking,
    inputs=[modelFolderVariable, ligandFolderVariable, dockingCenterVariable, innerBoxVariable],
    variables=[],
    outputs=[outputJobsVariable],
)
