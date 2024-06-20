import os
from flask import jsonify, request  # type: ignore
from HorusAPI import PluginPage, PluginEndpoint
import typing as t

docking_page = PluginPage(
    id="docking_analysis",
    name="Docking analysis",
    description="Analyze Glide results",
    html="docking_analysis.html",
)


def get_models():

    try:
        import prepare_proteins
        import pandas as pd

        import scipy

        import scipy.cluster

        data = request.get_json()

        docking_folder: t.Union[None, str] = data["dockingFolder"]
        ligand_separator: t.Union[None, str] = data["ligandSeparator"]

        if not ligand_separator:
            ligand_separator = "-"

        if docking_folder is None or not os.path.isdir(docking_folder):
            raise Exception("No valid docking folder selected")

        # The models are inside a folder called "input_models" inside the docking folder
        models_folder = os.path.join(docking_folder, "input_models")

        if models_folder is None or not os.path.isdir(models_folder):
            raise Exception("No 'input_models' folder found in the selected docking folder")

        # Loop over the input models folder to find the ligands
        output_models_folder = os.path.join(docking_folder, "output_models")

        if output_models_folder is None or not os.path.isdir(output_models_folder):
            raise Exception("No 'output_models' folder found in the selected docking folder")

        docking_results_dict = {}
        for imodel in os.listdir(output_models_folder):
            full_path = os.path.join(output_models_folder, imodel)

            ligands_for_model = {}
            # List the folder in search for every ligand (find the .csv for all of them)
            for f in os.listdir(full_path):
                if f.endswith(".csv") and not f.endswith("skip.csv"):
                    results_csv = ""
                    current_ligand_dict = {}

                    file_name, extension = os.path.splitext(f)
                    ligand_name = file_name.split(ligand_separator)[-1]
                    results_csv = os.path.join(
                        full_path, f"{imodel}{ligand_separator}{ligand_name}.csv"
                    )
                    if not os.path.isfile(results_csv):
                        raise ValueError(
                            f"Could not get file name using the provided separator. Trying to split '{file_name}' with '{ligand_separator}'"
                        )

                    # Save the info
                    current_ligand_dict["ligandName"] = ligand_name
                    current_ligand_dict["resultsCSVPath"] = results_csv
                    current_ligand_dict["logPath"] = os.path.join(
                        full_path, f"{imodel}{ligand_separator}{ligand_name}.log"
                    )
                    current_ligand_dict["logSubjobsPath"] = os.path.join(
                        full_path, f"{imodel}{ligand_separator}{ligand_name}_subjobs.log"
                    )

                    # Read the csv
                    current_model_csv = pd.read_csv(results_csv)

                    # Add an ID to each ligand
                    # Assign a unique id to each entry
                    current_model_csv["id"] = current_model_csv.index

                    current_ligand_dict["csvData"] = {
                        "rows": current_model_csv.to_dict(orient="records"),
                        "columns": current_model_csv.columns.tolist(),
                    }

                    ligands_for_model[ligand_name] = current_ligand_dict

            docking_results_dict[imodel] = ligands_for_model

        # models = prepare_proteins.proteinModels(models_folder)

        return jsonify({"ok": True, "results": docking_results_dict})
    except Exception as exc:
        return jsonify({"ok": False, "msg": str(exc)})


# Add the endpoint to the PluginPage
get_models_endpoint = PluginEndpoint(
    url="/api/get_models",
    methods=["POST"],
    function=get_models,
)

# Add the endpoint to the page
docking_page.addEndpoint(get_models_endpoint)


def extract_all_poses():

    import prepare_proteins

    data = request.get_json()
    docking_folder = data["dockingFolder"]
    ligand_name = data["ligandName"]

    # The models are inside a folder called "input_models" inside the docking folder
    models_folder = os.path.join(docking_folder, "input_models")

    models = prepare_proteins.proteinModels(models_folder)

    models.extractDockingPoses()


# Add the endpoint to the PluginPage
get_file_contents_endpoint = PluginEndpoint(
    url="/api/extract_all",
    methods=["POST"],
    function=extract_all_poses,
)

# Add the endpoint to the page
docking_page.addEndpoint(get_file_contents_endpoint)


def get_file_contents():
    try:
        file_path = request.args.get("path", None)

        if file_path is None:
            raise ValueError("No path provided")

        with open(file_path, "r") as f:
            read_file = f.read()
            return read_file
    except Exception as exc:
        return str(exc)


# Add the endpoint to the PluginPage
get_file_contents_endpoint = PluginEndpoint(
    url="/api/get_file_contents",
    methods=["GET"],
    function=get_file_contents,
)

# Add the endpoint to the page
docking_page.addEndpoint(get_file_contents_endpoint)


def analyze_atom_pairs():
    import prepare_proteins

    data = request.get_json()

    models_folder = data["modelsFolder"]
    docking_folder = data["dockingFolder"]
    atom_pairs_selection = data["atomPairs"]

    if models_folder is None or not os.path.isdir(models_folder):
        raise Exception("No valid model folder selected")

    if docking_folder is None or not os.path.isdir(docking_folder):
        raise Exception("No valid docking folder selected")

    if atom_pairs_selection is None or len(atom_pairs_selection) < 1:
        raise Exception("No atom pairs selected")

    models = prepare_proteins.proteinModels(models_folder)

    print(models.models_names)
    print(models.docking_ligands)

    atom_pairs = {}
    atom_pairs_for_pele = {}
    groups = []
    for model in models:
        atom_pairs[model] = {}
        atom_pairs_for_pele[model] = {}
        for selection in atom_pairs_selection:
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

    print(f"Atom pairs: {atom_pairs}")

    if atom_pairs == {}:
        raise Exception("No atom pairs were given, check the configuration of the block.")

    return jsonify({"dockingFolder": docking_folder, "modelsFolder": models_folder})


analyze_atom_pairs_endpoint = PluginEndpoint(
    url="/api/analyse_atoms",
    methods=["POST"],
    function=analyze_atom_pairs,
)

docking_page.addEndpoint(analyze_atom_pairs_endpoint)
