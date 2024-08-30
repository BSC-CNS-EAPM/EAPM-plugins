"""
Module containing the rDock to PELE block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
dockedFolder = PluginVariable(
    name="Docked folder",
    id="docked_folder",
    description="Folder where the docked ligands files have been stored.",
    type=VariableTypes.FOLDER,
)

# ==========================#
# Variable outputs
# ==========================#

modelsFolder = PluginVariable(
    name="Models folder",
    id="models_folder",
    description="Folder where the models will be stored.",
    type=VariableTypes.FOLDER,
)

##############################
#       Other variables      #
##############################

modelsToKeep = PluginVariable(
    name="Models to keep",
    id="models_to_keep",
    description="Number of the docked models with lowest scores to keep.",
    type=VariableTypes.INTEGER,
    defaultValue=1,
)


# Align action block
def rbdockToPELE(block: PluginBlock):

    import os
    import csv
    import pandas as pd
    import shutil
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import os
    from Bio import PDB

    def _rDockToDataFrame(folder_path):
        """
        Extracts the rDock data from the files in the folder and saves it in a DataFrame.

        Parameters:
        -----------
        folder_path : str
            The path to the folder containing the rDock files.
        """

        data = []
        storage_path = os.path.dirname(folder_path)

        for filename in [x for x in os.listdir(docked_folder) if x.endswith((".sd", ".sdf"))]:
            file_path = os.path.join(folder_path, filename)

            counter = 1
            score_bool = False
            conformer_bool = False

            # Open the file
            with open(file_path, "r") as file:
                for line in file:
                    if score_bool:
                        score = line.split()[0]
                    if conformer_bool:
                        ligand, conformer = line.split("-")
                        data.append([filename, counter, ligand, conformer, score])
                    if "$$$$" in line:
                        counter += 1
                    if ">  <SCORE>" in line:
                        score_bool = True
                    else:
                        score_bool = False
                    if ">  <s_lp_Variant>" in line:
                        conformer_bool = True
                    else:
                        conformer_bool = False

        # Write the extracted data to a CSV file
        output_file = "rDock_data.csv"
        with open(os.path.join(storage_path, output_file), "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["file_name", "file_entry", "ligand", "conformer", "rdock_score"])
            writer.writerows(data)

        print(" - rDock data extraction completed.")
        print(" - Data saved in {}".format(os.path.join(storage_path, output_file)))

    def _sdfSplitter(docked_folder_path):
        """
        If not done previously, split the rDock's output sdf files into its components.
        """

        path_docking = docked_folder_path
        path_docked_ligands = os.path.join(path_docking, "ligand_selection")

        # Generating storage folders
        if not os.path.isdir(path_docked_ligands):
            os.mkdir(path_docked_ligands)

        print(" - Splitting all the outputted sdf files.")

        with open(os.path.join(path_docking, ".docking.info"), "r") as f:
            for line in f:
                n_dockings = int(line.split()[10])

        # sdf splitting
        for file in [x for x in os.listdir(path_docking) if x.endswith((".sd", ".sdf"))]:
            input_file = os.path.join(path_docking, file)

            with open(input_file, "r") as file:
                sdf_content = file.read()

            previous_entry_id = None
            entries = sdf_content.strip().split("$$$$")
            auxiliar_counter = 0

            for entry in entries[:-1]:
                auxiliar_counter += 1
                count = auxiliar_counter - n_dockings * ((auxiliar_counter - 1) // n_dockings)
                entry = entry.strip()

                if entry:
                    entry_lines = entry.split("\n")
                    entry_id = None
                    for i in range(len(entry_lines)):
                        if entry_lines[i].strip().startswith(">  <s_lp_Variant>"):
                            entry_id = entry_lines[i + 1].strip()
                            break

                    if entry_id == previous_entry_id:
                        previous_entry_id = entry_id
                        output_file = os.path.join(
                            path_docked_ligands,
                            "{entry_id}_{c}.sdf".format(entry_id=entry_id, c=count),
                        )

                        with open(output_file, "w") as outfile:
                            outfile.write(entry)
                    else:
                        previous_entry_id = entry_id
                        output_file = os.path.join(
                            path_docked_ligands,
                            "{entry_id}_{c}.sdf".format(entry_id=entry_id, c=count),
                        )
                        with open(output_file, "w") as outfile:
                            outfile.write(entry)

    def _rDockPoseRetriever(folder_path, number_of_models):
        """
        Function to retrieve the ligands alone from the docked files
        """

        path = os.path.dirname(folder_path)
        df = pd.read_csv(os.path.join(path, "rDock_data.csv"))
        list_of_files = []
        for ligand in sorted(df.ligand.unique()):
            if number_of_models > len(df.ligand):
                raise Exception(
                    "The number of models to keep is higher than the number of models available."
                )

            df_indexed = df[df.ligand == ligand].reset_index(drop=True)
            df_trimmed = df_indexed.sort_values(by="rdock_score").head(number_of_models)

            df_trimmed["file_name"] = df_trimmed.apply(
                lambda row: f"{row['ligand']}-{row['conformer']}_{int(row.name) + 1}.sdf",
                axis=1,
            )

            list_of_files.extend(df_trimmed.file_name.tolist())

        return list_of_files

    def _rDocktoPELEModels(folder_path, models_folder, list_of_ligands):

        def _PDBConversor(file_in, path_out):
            """
            Converts the input file to PDB format using RDKit.

            Parameters
            ==========
            file_in : str
                Name of the file to be converted.
            path_out : str
                Path to store the output file.
            """

            file_name, file_format = file_in.split(".")
            file_name = os.path.basename(file_name)

            # Load the molecule depending on the input format
            mol = None
            if file_format == "sdf":
                suppl = Chem.SDMolSupplier(file_in)
                mol = suppl[0] if suppl else None
            elif file_format == "mol2":
                mol = Chem.MolFromMol2File(file_in)
            else:
                raise ValueError("Unsupported file format: {}".format(file_format))

            if mol is None:
                raise ValueError("Could not load molecule from file: {}".format(file_in))

            # Generate 3D coordinates if they don't exist
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())

            # Write the molecule to a PDB file
            output_path = os.path.join(path_out, "{}.pdb".format(file_name))
            Chem.MolToPDBFile(mol, output_path)

            return output_path

        def _PDBMerger(receptor, ligand):
            # Save the combined trajectory to a new PDB file
            parser = PDB.PDBParser(QUIET=True)

            # Parse both PDB files
            structure1 = parser.get_structure("structure1", receptor)
            structure2 = parser.get_structure("structure2", ligand)

            # Create a new structure object
            combined_structure = PDB.Structure.Structure("combined_structure")

            # Add models and chains from both structures to the new structure
            model_id = 0
            for model in structure1:
                model.id = model_id  # Ensure unique model ID
                combined_structure.add(model)
                model_id += 1

            for model in structure2:
                model.id = model_id  # Ensure unique model ID
                combined_structure.add(model)
                model_id += 1

            # Remove the original files
            os.remove(receptor)
            os.remove(ligand)

            # Write the combined structure to the ligand PDB file
            io = PDB.PDBIO()
            io.set_structure(combined_structure)
            io.save(ligand)

        # Getting paths
        docked_ligands_path = os.path.join(folder_path, "ligand_selection")
        with open(os.path.join(folder_path, ".docking.info"), "r") as f:
            for line in f:
                parameter_file = line.split()[6]
        with open(parameter_file) as file:
            for line in file:
                if "RECEPTOR_FILE" in line:
                    receptor_path = line.split()[1]
                    break

        path_to_receptor = os.path.dirname(receptor_path)

        # Generating PELE models
        print(" - Converting receptor from mol2 to pdb.")
        output_receptor_path = _PDBConversor(receptor_path, path_to_receptor)
        converted_receptor = os.path.basename(output_receptor_path)

        # Convert all the ligands to pdb and store them
        for ligand in os.listdir(docked_ligands_path):
            if ligand in list_of_ligands:
                ligand_name = ligand.split(".")[0]
                tmp_ligand_folder = os.path.join(models_folder, ligand_name)

                os.makedirs(tmp_ligand_folder, exist_ok=True)
                _ = _PDBConversor(
                    os.path.join(docked_ligands_path, ligand),
                    tmp_ligand_folder,
                )
                shutil.copy(
                    output_receptor_path,
                    tmp_ligand_folder,
                )

    if block.remote.name != "Local":
        raise Exception("This block is only available for local execution.")

    # Loading plugin variables
    docked_folder = block.inputs.get(dockedFolder.id, None)
    if docked_folder is None:
        raise Exception("The docked folder is not defined.")
    if not os.path.exists(docked_folder):
        raise Exception(f"The docked folder '{docked_folder}' does not exist.")
    if len([x for x in os.listdir(docked_folder) if x.endswith((".sd", ".sdf"))]) == 0:
        raise Exception(f"The docked folder '{docked_folder}' is empty.")

    models_to_keep = block.variables.get(modelsToKeep.id, 1)
    if not 0 < models_to_keep:
        raise Exception("The models to keep must be more than 0.")

    # Generatig models folder
    path = os.path.dirname(docked_folder)
    models_folder = os.path.join(path, "models")
    os.makedirs(models_folder, exist_ok=True)

    # Extracting rDock data
    _rDockToDataFrame(docked_folder)
    _sdfSplitter(docked_folder)
    list_of_files = _rDockPoseRetriever(docked_folder, block.variables.get(modelsToKeep.id, 1))
    _rDocktoPELEModels(docked_folder, models_folder, list_of_files)


rbDockToPELEBlock = PluginBlock(
    name="rDock to PELE",
    id="rbdock_to_pele",
    description="Extract information from rDock docking files to generate PELE models.",
    action=rbdockToPELE,
    variables=[modelsToKeep],
    inputs=[dockedFolder],
    outputs=[modelsFolder],
)
