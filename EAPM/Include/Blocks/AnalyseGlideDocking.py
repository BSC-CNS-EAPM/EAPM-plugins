import os
import datetime

from HorusAPI import (
    PluginBlock,
    VariableTypes,
    PluginVariable,
    VariableList,
    VariableGroup,
    Extensions,
)

# Input variables
dockingFolderVariable = PluginVariable(
    id="docking_folder",
    name="Docking Folder",
    description="Folder containing the docking results",
    type=VariableTypes.FOLDER,
    defaultValue="docking",
)

modelFolderVariable = PluginVariable(
    id="model_folder",
    name="Model folder",
    description="Folder containing the models",
    type=VariableTypes.FOLDER,
    defaultValue="models",
)

glideOutputVariable = PluginVariable(
    id="glide_output",
    name="Glide output",
    description="Glide output from the BSC calculations block",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_results"],
)

folderVariableGroup = VariableGroup(
    id="folder_variable_group",
    name="Folder variable group",
    description="Input the model and docking folders after a Docking Grid setup has been run",
    variables=[modelFolderVariable, dockingFolderVariable],
)

glideOutputVariableGroup = VariableGroup(
    id="glide_output_variable_group",
    name="Glide output variable group",
    description="Input the Glide output from the BSC calculations block",
    variables=[glideOutputVariable],
)

# Variables
ligandSeparatorVariable = PluginVariable(
    id="ligand_separator",
    name="Ligand separator",
    description="Separator used to separate the ligand name from the docking pose",
    type=VariableTypes.STRING,
    defaultValue="-",
)

atom1SelectionVariable = PluginVariable(
    id="protein_atom",
    name="Protein atom",
    description="Select an atom on the protein to calculate the distance",
    type=VariableTypes.ATOM,
)

atom2SelectionVariable = PluginVariable(
    id="ligand_atom",
    name="Ligand atom",
    description="Select an atom on the ligand to calculate the distance",
    type=VariableTypes.ATOM,
)

atomSelectionsGroupVariable = PluginVariable(
    id="group",
    name="Group",
    description="Name of the group selection to clusterize",
    type=VariableTypes.STRING,
)

overrideLigandNameVariable = PluginVariable(
    id="override_ligand_name",
    name="Override ligand name",
    description="Override the ligand name with a custom name",
    type=VariableTypes.STRING,
)

selectionsListVariable = VariableList(
    id="selections_list",
    name="Selections",
    description="List of selections to analyse",
    prototypes=[
        atom1SelectionVariable,
        atom2SelectionVariable,
        atomSelectionsGroupVariable,
        overrideLigandNameVariable,
    ],
)

posesFolderNameVariable = PluginVariable(
    id="poses_folder_name",
    name="Poses folder name",
    description="Name of the folder containing the docking poses",
    type=VariableTypes.STRING,
    defaultValue="best_docking_poses",
)

maxThresholdVariable = PluginVariable(
    id="max_threshold",
    name="Max threshold",
    description="Maximum threshold to consider a pose as a good pose",
    type=VariableTypes.FLOAT,
    defaultValue=5,
)

# Output variables
outputPosesVariable = PluginVariable(
    id="output_poses",
    name="Output poses",
    description="Output folder containing the best docking poses",
    type=VariableTypes.FOLDER,
)

analyseGlideOutputVariable = PluginVariable(
    id="glide_results_output",
    name="Glide results output",
    description="Output results of the Glide analysis",
    type=VariableTypes.CUSTOM,
    allowedValues=["glide_output"],
)


def analyseDockingAction(block: PluginBlock):
    if block.selectedInputGroup == "folder_variable_group":
        folder_to_analyse = block.inputs.get("docking_folder", "docking")
        model_folder = block.inputs.get("model_folder", "models")
    else:
        bsc_result = block.inputs.get("glide_output")
        folder_to_analyse = bsc_result["dock_folder"]
        model_folder = bsc_result["model_folder"]

    if folder_to_analyse is None or not os.path.isdir(folder_to_analyse):
        raise Exception("No valid docking folder selected")

    if model_folder is None or not os.path.isdir(model_folder):
        raise Exception("No valid model folder selected")

    import prepare_proteins

    # Get the docking results
    models = prepare_proteins.proteinModels(model_folder)

    # Generate the atom pairs based on the selections
    selections = block.variables.get("selections_list", [])
    output_poses = block.variables.get("poses_folder_name", "best_docking_poses")

    # If the folder to analyse is not on the current working directory, we need to
    # copy it to the curret working directory
    current_folder = os.getcwd()

    # If the folder starts with a different path, we need to copy it to the current dir
    # If its a relative path (under the current dir) pass
    if (
        not folder_to_analyse.startswith(current_folder)
        and not os.path.basename(folder_to_analyse) == folder_to_analyse
    ):
        # Copy the folder to the current working directory
        import shutil

        shutil.copytree(
            folder_to_analyse, os.path.join(current_folder, folder_to_analyse), dirs_exist_ok=True
        )

    folder_to_analyse = os.path.abspath(folder_to_analyse)

    # If we are on the remote, create the output folder
    remote_folder = str(block.flow.savedID + "_" + str(datetime.datetime.now().timestamp()))
    remote_folder = os.path.join(block.remote.workDir, remote_folder)

    block.extraData["remote_folder"] = remote_folder

    # Create the remote folder
    block.remote.remoteCommand(f"mkdir -p {remote_folder}")

    separatorValue = block.variables.get("ligand_separator", "-")

    if selections is None or selections == []:

        # Extract single docking poses
        docking_data = models.getSingleDockingData(
            models.models_names[0], models.docking_ligands[models.models_names[0]]
        )
        models.extractDockingPoses(
            docking_data, folder_to_analyse, output_folder="poses", separator=separatorValue
        )

        extractDockingPoses(
            block,
            models,
            docking_data,
            block.extraData["final_path_folder_to_analyse"],
            output_poses,
        )

        block.setOutput("output_poses", output_poses)

        glideOutput = {
            "poses_folder": output_poses,
            "models_folder": model_folder,
            "atom_pairs": atom_pairs,
        }

        block.setOutput("glide_results_output", glideOutput)

        print("Successfully extracted a single docking pose")

        return

    atom_pairs = {}
    atom_pairs_for_pele = {}
    groups = []
    for model in models:
        atom_pairs[model] = {}
        atom_pairs_for_pele[model] = {}
        for selection in selections:
            current_group = selection["group"]
            if current_group not in groups:
                groups.append(current_group)
            atom1 = selection["protein_atom"]

            protein_chain = atom1["chainID"]
            protein_resnum = atom1["residue"]
            protein_atom = atom1["auth_atom_id"]

            protein_tuple = (protein_chain, protein_resnum, protein_atom)

            atom2 = selection["ligand_atom"]
            ligandName = atom2["auth_comp_id"]

            ligand_chain = atom2["chainID"]
            ligand_resnum = atom2["residue"]
            ligand_atom = atom2["auth_atom_id"]

            atom_tuple = (ligand_chain, ligand_resnum, ligand_atom)

            if (
                selection.get("override_ligand_name") is not None
                and selection["override_ligand_name"] != ""
            ):
                ligandName = selection["override_ligand_name"]

            if atom_pairs[model].get(ligandName, None) is None:
                atom_pairs[model][ligandName] = []
                atom_pairs_for_pele[model][ligandName] = []

            atom_pairs[model][ligandName].append((protein_tuple, ligand_atom))
            atom_pairs_for_pele[model][ligandName].append((protein_tuple, atom_tuple))

    if atom_pairs == {}:
        raise Exception("No atom pairs were given, check the configuration of the block.")

    print("Starting docking analysis")

    analyseDocking(
        block,
        models,
        folder_to_analyse,
        atom_pairs=atom_pairs,
        return_failed=True,
        separator=separatorValue,
    )

    print("Docking analysis finished")

    metric_distances = {}  # Define the global dictionary
    for group in groups:
        metric_distances[group] = {}  # Define the metric nested dictionary
        for model in models:
            metric_distances[group][model] = {}  # Define the model nested dictionary
            for ligand in models.docking_ligands[model]:
                # Define the ligand nested dictionary with all the docking distances list
                metric_distances[group][model][ligand] = models.getDockingDistances(model, ligand)

    models.combineDockingDistancesIntoMetrics(metric_distances, overwrite=True)

    # Generate an HTML report containing the results
    if models.docking_data is not None:
        html = models.docking_data.to_html()

        Extensions().storeExtensionResults(
            "horus", "html_loader", data={"html": html}, title="Docking results"
        )

    max_threshold = float(block.variables.get("max_threshold", 5))

    best_poses = models.getBestDockingPosesIteratively(
        metric_distances, max_threshold=max_threshold
    )

    if len(best_poses) == 0:
        raise Exception("No best poses found with the given threshold. Try a lower threshold")

    html = best_poses.to_html()

    Extensions().storeExtensionResults(
        "horus", "html_loader", data={"html": html}, title="Best docking poses"
    )

    extractDockingPoses(
        block, models, best_poses, block.extraData["final_path_folder_to_analyse"], output_poses
    )

    print("Docking analysis finished")

    block.setOutput("output_poses", output_poses)

    glideOutput = {
        "poses_folder": output_poses,
        "models_folder": model_folder,
        "atom_pairs": atom_pairs,
    }

    block.setOutput("glide_results_output", glideOutput)


analyseGlideDocking = PluginBlock(
    name="Analyse Glide Docking",
    description="Analyse the docking results from Glide",
    inputGroups=[folderVariableGroup, glideOutputVariableGroup],
    variables=[
        ligandSeparatorVariable,
        maxThresholdVariable,
        posesFolderNameVariable,
        selectionsListVariable,
    ],
    outputs=[outputPosesVariable, analyseGlideOutputVariable],
    action=analyseDockingAction,
)


def analyseDocking(
    block: PluginBlock,
    models,
    docking_folder,
    protein_atoms=None,
    atom_pairs=None,
    skip_chains=False,
    return_failed=False,
    ignore_hydrogens=False,
    separator="-",
    overwrite=True,
):
    """
    Analyse a Glide Docking simulation. The function allows to calculate ligand
    distances with the options protein_atoms or protein_pairs. With the first option
    the analysis will calculate the closest distance between the protein atoms given
    and any ligand atom (or heavy atom if ignore_hydrogens=True). The analysis will
    also return which ligand atom is the closest for each pose. On the other hand, with
    the atom_pairs option only distances for the specific atom pairs between the
    protein and the ligand will be calculated.

    The protein_atoms dictionary must contain as keys the model names (see iterable of this class),
    and as values a list of tuples, with each tuple representing a protein atom:
        {model1_name: [(chain1_id, residue1_id, atom1_name), (chain2_id, residue2_id, atom2_name), ...], model2_name:...}

    The atom pairs must be given in a dictionary with each key representing the name
    of a model and each value  a sub dictionary with the ligands as keys and a list of the atom pairs
    to calculate in the format:
        {model1_name: { ligand_name : [((chain1_id, residue1_id, atom1_name), (chain2_id, residue2_id, atom2_name)), ...],...} model2_name:...}

    Parameters
    ===========
    docking_folder : str
        Path to the folder where the docking results are (the format comes from the setUpGlideDocking() function.
    protein_atoms : dict
        Protein atoms to use for the closest distance calculation.
    atom_pairs : dict
        Protein and ligand atoms to use for distances calculation.
    skip_chains : bool
        Consider chains when atom tuples are given?
    return_failed : bool
        Return failed dockings as a list?
    ignore_hydrogens : bool
        With this option ligand hydrogens will be ignored for the closest distance (i.e., protein_atoms) calculation.
    separator : str
        Symbol to use for separating protein from ligand names. Should not be found in any model or ligand name.
    overwrite : bool
        Rerun analysis.
    """

    import prepare_proteins
    import json
    import pandas as pd

    # Create analysis folder
    if not os.path.exists(docking_folder + "/.analysis"):
        os.mkdir(docking_folder + "/.analysis")

    # Create analysis folder
    if not os.path.exists(docking_folder + "/.analysis/atom_pairs"):
        os.mkdir(docking_folder + "/.analysis/atom_pairs")

    # Copy analyse docking script (it depends on Schrodinger Python API so we leave it out to minimise dependencies)
    prepare_proteins._copyScriptFile(docking_folder + "/.analysis", "analyse_docking.py")
    script_path = docking_folder + "/.analysis/._analyse_docking.py"

    # Write protein_atoms dictionary to json file
    if protein_atoms != None:
        with open(docking_folder + "/.analysis/._protein_atoms.json", "w") as jf:
            json.dump(protein_atoms, jf)

    # Write atom_pairs dictionary to json file
    if atom_pairs != None:
        with open(docking_folder + "/.analysis/._atom_pairs.json", "w") as jf:
            json.dump(atom_pairs, jf)

    # Upload the folder to the remote
    remote_docking_folder = block.remote.sendData(
        docking_folder, block.extraData["remote_folder"]
    )

    block.extraData["final_path_folder_to_analyse"] = remote_docking_folder

    # Replace the os.system call with a remote call
    # os.system(command)

    # Test if we have valid schrodinger installation
    schrodinger = block.remote.remoteCommand("echo $SCHRODINGER")

    if schrodinger is None or schrodinger == "":
        raise Exception(f"No valid Schrodinger installation found on remote {block.remote.name}")
    else:
        print(f"Schrodinger installation found on remote {block.remote.name}: {schrodinger}")

    command = (
        f"{schrodinger}/run "
        + remote_docking_folder
        + "/.analysis/._analyse_docking.py "
        + remote_docking_folder
    )
    if atom_pairs != None:
        command += " --atom_pairs " + remote_docking_folder + "/.analysis/._atom_pairs.json"
    elif protein_atoms != None:
        command += " --protein_atoms " + remote_docking_folder + "/.analysis/._protein_atoms.json"
    if skip_chains:
        command += " --skip_chains"
    if return_failed:
        command += " --return_failed"
    if ignore_hydrogens:
        command += " --ignore_hydrogens"
    command += " --separator " + separator
    command += " --only_models " + ",".join(models.models_names)
    if overwrite:
        command += " --overwrite "

    # Execute the command
    output = block.remote.remoteCommand(command)
    print(output)

    # Download the results
    docking_folder = block.remote.getData(remote_docking_folder, os.getcwd())

    # Read the CSV file into pandas
    if not os.path.exists(docking_folder + "/.analysis/docking_data.csv"):
        raise ValueError(
            "Docking analysis failed. Check the output of the analyse_docking.py script."
        )

    models.docking_data = pd.read_csv(docking_folder + "/.analysis/docking_data.csv")
    # Create multiindex dataframe
    models.docking_data.set_index(["Protein", "Ligand", "Pose"], inplace=True)

    for f in os.listdir(docking_folder + "/.analysis/atom_pairs"):
        model = f.split(separator)[0]
        ligand = f.split(separator)[1].split(".")[0]

        # # Read the CSV file into pandas
        models.docking_distances.setdefault(model, {})
        models.docking_distances[model][ligand] = pd.read_csv(
            docking_folder + "/.analysis/atom_pairs/" + f
        )
        models.docking_distances[model][ligand].set_index(
            ["Protein", "Ligand", "Pose"], inplace=True
        )

        models.docking_ligands.setdefault(model, [])
        if ligand not in models.docking_ligands[model]:
            models.docking_ligands[model].append(ligand)

    if return_failed:
        with open(docking_folder + "/.analysis/._failed_dockings.json") as jifd:
            failed_dockings = json.load(jifd)
        return failed_dockings


def extractDockingPoses(
    block: PluginBlock,
    models,
    docking_data,
    docking_folder,
    output_folder,
    separator="-",
    only_extract_new=True,
    covalent_check=True,
    remove_previous=False,
):
    """
    Extract docking poses present in a docking_data dataframe. The docking DataFrame
    contains the same structure as the self.docking_data dataframe, parameter of
    this class. This dataframe makes reference to the docking_folder where the
    docking results are contained.

    Parameters
    ==========
    dockign_data : pandas.DataFrame
        Datframe containing the poses to be extracted
    docking_folder : str
        Path the folder containing the docking results
    output_folder : str
        Path to the folder where the docking structures will be saved.
    separator : str
        Symbol used to separate protein, ligand, and docking pose index.
    only_extract_new : bool
        Only extract models not present in the output_folder
    remove_previous : bool
        Remove all content in the output folder
    """

    import shutil

    # Check the separator is not in model or ligand names
    for model in models.docking_ligands:
        if separator in str(model):
            raise ValueError(
                "The separator %s was found in model name %s. Please use a different separator symbol."
                % (separator, model)
            )
        for ligand in models.docking_ligands[model]:
            if separator in ligand:
                raise ValueError(
                    "The separator %s was found in ligand name %s. Please use a different separator symbol."
                    % (separator, ligand)
                )

    # Remove output_folder
    if os.path.exists(output_folder):
        if remove_previous:
            shutil.rmtree(output_folder)

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    else:
        # Gather already extracted models
        if only_extract_new:
            extracted_models = set()
            for model in os.listdir(output_folder):
                if not os.path.isdir(output_folder + "/" + model):
                    continue
                for f in os.listdir(output_folder + "/" + model):
                    if f.endswith(".pdb"):
                        m, l = f.split(separator)[:2]
                        extracted_models.add((m, l))

            # Filter docking data to not include the already extracted models
            extracted_indexes = []
            for i in docking_data.index:
                if i[:2] in extracted_models:
                    extracted_indexes.append(i)
            docking_data = docking_data[~docking_data.index.isin(extracted_indexes)]
            if docking_data.empty:
                print("All models were already extracted!")
                print("Set only_extract_new=False to extract them again!")
                return
            else:
                print(f"{len(extracted_models)} models were already extracted!")
                print(f"Extracting {docking_data.shape[0]} new models")

    # Copy analyse docking script (it depends on schrodinger so we leave it out.)
    _copyScriptFile(output_folder, "extract_docking.py")
    script_path = output_folder + "/._extract_docking.py"

    output_folder_abs = os.path.abspath(output_folder)

    # Move to output folder
    os.chdir(output_folder)

    # Save given docking data to csv
    dd = docking_data.reset_index()
    dd.to_csv("._docking_data.csv", index=False)

    # Upload the data to the remote
    remote_folder = block.extraData["remote_folder"]

    # Create the remote folder
    block.remote.remoteCommand(f"mkdir -p {remote_folder}")

    # Upload the folder to the remote
    final_path = block.remote.sendData(output_folder_abs, remote_folder)

    print(f"Uploaded {output_folder} to {final_path}")
    schrodinger = block.remote.remoteCommand("echo $SCHRODINGER")

    # Execute docking analysis
    command = (
        f"cd {final_path} && "
        + f"{schrodinger}/run ._extract_docking.py ._docking_data.csv "
        + docking_folder
        + " --separator "
        + separator
    )

    # Execute the command
    output = block.remote.remoteCommand(command)
    print(output)

    # Download the results
    block.remote.getData(final_path, os.path.join(os.getcwd(), ".."))

    # Remove the remote folder
    output = block.remote.remoteCommand(f"rm -rf {remote_folder}")

    # Remove docking data
    os.remove("._docking_data.csv")

    # move back to folder
    os.chdir("..")

    # Check models for covalent residues
    for protein in os.listdir(output_folder):
        if not os.path.isdir(output_folder + "/" + protein):
            continue
        for f in os.listdir(output_folder + "/" + protein):
            if covalent_check:
                models._checkCovalentLigands(
                    protein, output_folder + "/" + protein + "/" + f, check_file=True
                )


def _copyScriptFile(output_folder, script_name, no_py=False, subfolder=None, hidden=True):
    """
    Copy a script file from the prepare_proteins package.

    Parameters
    ==========

    """
    from pkg_resources import resource_stream, Requirement
    import io

    # Get script
    path = "prepare_proteins/scripts"
    if subfolder != None:
        path = path + "/" + subfolder

    script_file = resource_stream(Requirement.parse("prepare_proteins"), path + "/" + script_name)
    script_file = io.TextIOWrapper(script_file)

    # Write control script to output folder
    if no_py == True:
        script_name = script_name.replace(".py", "")

    if hidden:
        output_path = output_folder + "/._" + script_name
    else:
        output_path = output_folder + "/" + script_name

    with open(output_path, "w") as sof:
        for l in script_file:
            sof.write(l)
