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
    import openbabel as ob
    import re

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
            output_path = os.path.join(path_out, "{}.pdb".format(file_name))

            if not file_format == "pdb":

                conv = ob.OBConversion()
                conv.SetInAndOutFormats(file_in.split(".")[-1], "pdb")
                mol = ob.OBMol()
                conv.ReadFile(mol, file_in)
                conv.WriteFile(mol, output_path)

            return output_path

        def _PDBMerger(receptor, ligand):
            """
            Merges receptor and ligand into a single file with the characteristics required by
            PELE to work properly.

            Parameters
            ==========
            receptor : str
                Path to the receptor file we want to merge.
            ligand : str
                Path to the ligand file we want to merge.
            """

            def _receptorModifier(receptor):
                """
                Changes residues numbers and counts number of atoms in the receptor.

                Parameters
                ==========
                receptor : str
                    Path to the receptor file we want to merge.

                Returns
                =======
                file_mod_prot : str
                    Path to the modified version of the receptor.
                ligand_cont_num : int
                    Number of atoms in the receptor.
                """

                file_mod_prot = receptor.split(".pdb")[0] + "_mod.pdb"

                with open(receptor, "r") as pdb_file:
                    lines = pdb_file.readlines()

                modified_lines = []
                residue_index = 1000
                previous_residue = None
                previous_chain_residue = None

                # Changeing residue number
                for line in lines:

                    # Change lines of the receptor
                    if line.startswith("ATOM"):
                        residue_letters = line[17:20].strip()
                        chain_residue = line[21:26].strip()

                        if (residue_letters != previous_residue) or (
                            chain_residue != previous_chain_residue
                        ):
                            residue_index += 1

                        line = line[:22] + str(residue_index) + line[26:]
                        modified_lines.append(line)
                        previous_residue = residue_letters
                        previous_chain_residue = chain_residue

                    elif line.startswith("TER"):
                        pass

                    # Change lines of waters/metals
                    elif line.startswith("HETATM"):
                        residue_letters = line[17:20].strip()
                        chain_residue = line[21:26].strip()

                        if (residue_letters != previous_residue) or (
                            chain_residue != previous_chain_residue
                        ):
                            residue_index += 1

                        line = line[:22] + str(residue_index) + line[26:]
                        modified_lines.append(line)
                        previous_residue = residue_letters
                        previous_chain_residue = chain_residue

                    else:
                        pass

                # Write the modified lines back to the PDB file
                with open(file_mod_prot, "w") as pdb_file:
                    pdb_file.writelines(modified_lines)

                # Counting number of atoms of the protein
                ligand_cont_num = 0

                with open(file_mod_prot) as filein:
                    for line in filein:
                        sline = line.split()
                        if sline[0] == "ATOM":
                            ligand_cont_num += 1

                return file_mod_prot, ligand_cont_num

            def _ligandAtomChainNumberModifier(ligand):
                """
                Modifies the ligand file to have the characteristics needed for the PELE simulation:
                like the ligand chain or the chain name.

                Parameters
                ==========
                ligand : str
                    Path to the ligand file we want to merge.

                Returns
                =======
                file_mod_lig : str
                    Path to the modified version of the ligand.
                """

                file_mod1_lig = ligand.split(".pdb")[0] + "_m.pdb"
                file_mod_lig = ligand.split(".pdb")[0] + "_mod.pdb"

                # Modifying the atom names
                cont = 1
                with open(ligand) as filein:
                    with open(file_mod1_lig, "w") as fileout:
                        for line in filein:
                            if line.startswith("HETATM"):
                                sline = line.split()
                                beginning_of_line = line[0:13]
                                end_of_line = line[17:]
                                atom = "".join([c for c in sline[2] if c.isalpha()])
                                new_line = (
                                    beginning_of_line + (atom + str(cont)).ljust(4) + end_of_line
                                )
                                fileout.writelines(new_line)
                                cont += 1

                # Change the chain name and number
                with open(file_mod1_lig, "r") as filein:
                    with open(file_mod_lig, "w") as fileout:
                        for line in filein:
                            if re.search("UN.......", line):
                                line = re.sub("UN.......", "LIG L 900", line)
                                fileout.writelines(line)
                            else:
                                fileout.writelines(line)

                return file_mod_lig

            def _receptorLigandMerger(receptor, ligand):
                """
                Merges the receptor and the ligand into a single file.

                Parameters
                ==========
                receptor : str
                    Path to the receptor file we want to merge.
                ligand : str
                    Path to the ligand file we want to merge.

                Returns
                =======
                writing_path : str
                    Path to the newly created pdb with both molecules.
                """

                output_file = "intermediate.pdb"
                path_files = os.path.dirname(ligand)
                writing_path = os.path.join(path_files, output_file)

                # Joining ligand and receptor and directing warnings to out.txt
                os.system(
                    "obabel {receptor} {ligand} -O {output} 2> out.txt".format(
                        receptor=receptor, ligand=ligand, output=writing_path
                    )
                )

                return writing_path

            def _inputAdapter(merged, output_file, ligand_cont_num):
                """
                Changes the atom numeration of the ligand according to the number
                of atoms of the receptor.

                Parameters
                ==========
                merged : str
                    Path to the merged file.
                output_file : str
                    Path of the definitive merged file.
                ligand_cont_num : int
                    Number of atoms of the receptor.
                """

                # Changing the atom numeration of the ligand.
                with open(merged, "r") as filein:
                    with open(output_file, "w") as fileout:
                        for line in filein:
                            sline = line.split()
                            if sline[0] == "ATOM" or sline[0] == "HETATM":
                                if re.search("HETATM.....", line):
                                    # Considering length of the protein
                                    if len(str(ligand_cont_num + 1)) < 5:
                                        line = re.sub(
                                            "HETATM.....",
                                            "HETATM " + str(ligand_cont_num + 1),
                                            line,
                                        )
                                    elif len(str(ligand_cont_num + 1)) == 5:
                                        line = re.sub(
                                            "HETATM.....",
                                            "HETATM" + str(ligand_cont_num + 1),
                                            line,
                                        )

                                    fileout.writelines(line)
                                    ligand_cont_num += 1
                                else:
                                    fileout.writelines(line)

            def _intermediateFilesRemover(output_file):
                """
                Removes all the intermediate files that have been
                generated in the process.

                Parameters
                ==========
                output_file : str
                    Path of the definitive merged file.
                """

                working_directory = os.path.dirname(output_file)
                input_PELE_file = os.path.basename(output_file)

                for file in os.listdir(working_directory):
                    if file != input_PELE_file:
                        os.remove(os.path.join(working_directory, file))

            ligand_name = os.path.basename(ligand)
            working_directory = os.path.dirname(ligand)
            output_file = os.path.join(working_directory, ligand_name)

            file_mod_prot, ligand_cont_num = _receptorModifier(receptor)
            file_mod_lig = _ligandAtomChainNumberModifier(ligand)
            intermediate = _receptorLigandMerger(file_mod_prot, file_mod_lig)
            _inputAdapter(intermediate, output_file, ligand_cont_num)
            _intermediateFilesRemover(output_file)

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
                tmp_ligand_path = os.path.join(tmp_ligand_folder, ligand.split(".")[0] + ".pdb")
                tmp_receptor_path = os.path.join(tmp_ligand_folder, converted_receptor)

                os.makedirs(tmp_ligand_folder, exist_ok=True)
                _ = _PDBConversor(
                    os.path.join(docked_ligands_path, ligand),
                    tmp_ligand_folder,
                )
                shutil.copy(
                    output_receptor_path,
                    tmp_ligand_folder,
                )

                _PDBMerger(tmp_receptor_path, tmp_ligand_path)

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
    if len(os.listdir(models_folder)) != 0:
        print("Warning: Deleting previous models folder")
        for model in os.listdir(models_folder):
            shutil.rmtree(os.path.join(models_folder, model))

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
