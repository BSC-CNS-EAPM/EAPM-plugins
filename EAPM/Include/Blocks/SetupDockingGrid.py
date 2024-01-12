import os

from HorusAPI import InputBlock, SlurmBlock, VariableTypes, PluginVariable, VariableGroup

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
gridOutputVariable = PluginVariable(
    id="grid_output",
    name="GRID",
    description="Grid to be used in the docking calculation",
    type=VariableTypes.CUSTOM,
    allowedValues=["grid_output"],
)


# Action
def glideDocking(block: SlurmBlock):
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

    if block.selectedInputGroup != "single_model":
        # Get the common residues
        common_residues = block.inputs.get("multimodel_common_residue", {})

        # Parse the atomcenter for each model
        center_atoms = {}  # Create dictionary to store the atom 3-element tuple for each model
        for model in models:  # Iterate the models inside the library
            for r in models.structures[
                model
            ].get_residues():  # Iterate the residues for each Bio.PDB.Structure object
                if (
                    r.id[1] == common_residues[model][0]
                ):  # Check that the residue matches the defined index
                    center_atoms[model] = (
                        r.get_parent().id,
                        r.id[1],
                        RESIDUE_DICTIONARY[r.resname],
                    )  # Store the corresponsing tuple.
    else:
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

    from utils import launchCalculationAction

    launchCalculationAction(block, jobs, "glide", ["grid"])


def downloadGridResults(block: SlurmBlock):
    from utils import downloadResultsAction

    downloadResultsAction(block)

    # Set as the output the same models and ligands folder as the input
    grid_output = {
        "model_folder": block.inputs["model_folder"],
        "ligand_folder": block.inputs["ligand_folder"],
    }

    block.setOutput("grid_output", grid_output)


from utils import BSC_JOB_VARIABLES

setupDockingGrid = SlurmBlock(
    name="Setup Docking Grid",
    description="Prepare the files for a GLIDE docking calculation",
    initialAction=glideDocking,
    finalAction=downloadGridResults,
    inputGroups=[
        VariableGroup(
            id="single_model",
            name="Single model",
            description="Specify the box center and size for a single model",
            variables=[
                modelFolderVariable,
                ligandFolderVariable,
                dockingCenterVariable,
                innerBoxVariable,
            ],
        ),
        VariableGroup(
            id="multiple_models",
            name="Multiple models common residue",
            description="Provide a common residue for all the models",
            variables=[
                modelFolderVariable,
                ligandFolderVariable,
                dockingCenterVariable,
                innerBoxVariable,
                PluginVariable(
                    id="multimodel_common_residue",
                    name="Common residue",
                    description="Select a residue that is common to all the models",
                    type=VariableTypes.CUSTOM,
                    allowedValues=["multimodel_common_residue_grid"],
                ),
            ],
        ),
    ],
    variables=BSC_JOB_VARIABLES,
    outputs=[gridOutputVariable],
)

RESIDUE_DICTIONARY = {
    "ALA": "CB",
    "CYS": "SG",
    "ASP": "CG",
    "GLU": "CD",
    "PHE": "CZ",
    "GLY": "CA",
    "HIS": "NE2",
    "ILE": "CB",
    "LYS": "NZ",
    "LEU": "CG",
    "MET": "SD",
    "ASN": "CG",
    "PRO": "CG",
    "GLN": "CD",
    "ARG": "CZ",
    "SER": "OG",
    "THR": "OG1",
    "TRP": "NE1",
    "VAL": "CB",
    "TYR": "OH",
}
