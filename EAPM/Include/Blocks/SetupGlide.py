import os
import shutil

from HorusAPI import SlurmBlock, VariableTypes, PluginVariable, VariableGroup

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

folderInputGroup = VariableGroup(
    id="folder_input_group",
    name="Folder input group",
    description="Input the model and ligand folders after a Dcoking Grid setup has been run",
    variables=[modelFolderVariable, ligandFolderVariable],
)

gridOutputVariable = PluginVariable(
    id="grid_output",
    name="Grid output",
    description="Grid calculation output from the BSC calculations block",
    type=VariableTypes.CUSTOM,
    allowedValues=["grid_output"],
)

# Other variables
posesPerLigandVariable = PluginVariable(
    id="poses_per_ligand",
    name="Poses per ligand",
    description="Number of poses to generate per ligand",
    type=VariableTypes.INTEGER,
    defaultValue=10000,
)

outputDockingResultsVariable = PluginVariable(
    id="docking_results",
    name="Docking results",
    description="Information containing the docking results for the analyse block",
    type=VariableTypes.CUSTOM,
    allowedValues=["glide_docking_results"],
)


def setupGlideDocking(block: SlurmBlock):
    import prepare_proteins

    if block.selectedInputGroup == "folder_input_group":
        models_folder = block.inputs.get("model_folder")
        ligand_folder = block.inputs.get("ligand_folder")
    else:
        grid_output = block.inputs.get("grid_output")

        if grid_output is None:
            raise Exception("No valid grid output selected")

        models_folder = grid_output.get("model_folder")
        ligand_folder = grid_output.get("ligand_folder")

        if models_folder is None or not os.path.isdir(models_folder):
            raise Exception("No valid input")

    if models_folder is None or not os.path.isdir(models_folder):
        raise Exception("No valid models folder selected")

    if ligand_folder is None or not os.path.isdir(ligand_folder):
        raise Exception("No valid ligands folder selected")

    # If no maes are inside the ligan folder, raise an exception
    mae_found = False
    for file in os.listdir(ligand_folder):
        if file.endswith(".mae"):
            mae_found = True
            break

    if not mae_found:
        raise Exception(f"No .mae files found in {ligand_folder}")

    # Get the original pdb models
    original_pdb_folder = os.path.basename(models_folder).replace("_mae", "")

    if not os.path.isdir(original_pdb_folder):
        raise Exception(
            f"PDB models folder ({original_pdb_folder}) not found. Please convert the models to PDB using the 'MAE to PDB' block and keep the original PDB models folder"
        )

    block.extraData["original_pdb_folder"] = original_pdb_folder
    block.extraData["ligand_folder"] = ligand_folder
    block.extraData["models_folder"] = original_pdb_folder

    models = prepare_proteins.proteinModels(original_pdb_folder)

    poses_per_lig = block.variables.get("poses_per_ligand", 10000)

    # Set the ligand folder path relative to the current directory
    relative_ligand_folder = os.path.basename(ligand_folder)

    # If the ligand folder is not in the current directory, copy it
    if ligand_folder != os.path.join(os.getcwd(), relative_ligand_folder):
        print(f"Copying ligand folder {ligand_folder} to flow directory")
        shutil.copytree(ligand_folder, os.path.join(os.getcwd(), ligand_folder))

    jobs = models.setUpGlideDocking(
        "docking", "grid", relative_ligand_folder, poses_per_lig=poses_per_lig
    )

    jobs_created = len(jobs)

    if jobs_created == 0:
        raise Exception("No jobs created. Did the Glide Grid block produce the correct output?")

    from utils import launchCalculationAction

    launchCalculationAction(
        block, jobs, "glide", uploadFolders=["docking", "grid", relative_ligand_folder]
    )


def downloadGlideDocking(block: SlurmBlock):
    from utils import downloadResultsAction

    downloadResultsAction(block)

    results_data = {
        "ligand_folder": block.extraData["ligand_folder"],
        "model_folder": block.extraData["models_folder"],
        "dock_folder": "docking",
    }

    # Check on the output logs for each model if there was an error
    for model in os.listdir(os.path.join(os.getcwd(), "docking", "output_models")):
        # Skip .DS_Store files
        if model.startswith("."):
            continue

        # Get the log file (the only file that ends with .log)
        log_file = None
        for file in os.listdir(os.path.join(os.getcwd(), "docking", "output_models", model)):
            if file.endswith(".log"):
                log_file = os.path.join(os.getcwd(), "docking", "output_models", model, file)
                break

        if log_file is None:
            continue

        if not os.path.isfile(log_file):
            continue

        with open(log_file, "r") as f:
            log_content = f.read()

            if "ROUGH POSE REFINE FAILED" in log_content:
                raise Exception(
                    f"Rough pose refine failed for model {model}, try a different grid size"
                )

    block.setOutput(outputDockingResultsVariable.id, results_data)


from utils import BSC_JOB_VARIABLES

block_variables = BSC_JOB_VARIABLES + [posesPerLigandVariable]

setupGlideBlock = SlurmBlock(
    name="Run Glide",
    description="Run a glide docking calculation.",
    inputGroups=[
        folderInputGroup,
        VariableGroup(
            id="docking_grid_group",
            name="Docking grid",
            description="Input the docking grid output from the BSC calculations block",
            variables=[gridOutputVariable],
        ),
    ],
    variables=block_variables,
    initialAction=setupGlideDocking,
    finalAction=downloadGlideDocking,
    outputs=[outputDockingResultsVariable],
)
