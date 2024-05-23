import json
import os
import shutil

import bsc_calculations
import pandas as pd
import prepare_proteins

from HorusAPI import PluginBlock, PluginVariable, VariableGroup, VariableTypes

# TODO  Configure the inputs correctly

# ==========================#
# Variable inputs
# ==========================#
glideOutputVariable = PluginVariable(
    id="glide_output",
    name="Glide output",
    description="Glide output from the BSC calculations block",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_results"],
)
conservedResidues = PluginVariable(
    name="Conserved residues",
    id="conserved_indexes",
    description="The conserved residues",
    type=VariableTypes.CUSTOM,
    defaultValue=None,
)
residueProtein = PluginVariable(
    name="Atom Protein", id="resi_id1", description="atom1", type=VariableTypes.ATOM
)
residueLigand = PluginVariable(
    name="Atom Ligand", id="resi_id2", description="atom2", type=VariableTypes.ATOM
)
resNameProt = PluginVariable(
    name="Protein residue name",
    id="res_name_prot",
    description="The protein residue name",
    type=VariableTypes.STRING,
    defaultValue="CYS",
)
resNameLig = PluginVariable(
    name="Ligand residue name",
    id="res_name_lig",
    description="The ligand residue name",
    type=VariableTypes.STRING,
    defaultValue="SG",
)
ligandName = PluginVariable(
    name="Ligand name",
    id="ligand_name",
    description="The ligand name",
    type=VariableTypes.STRING,
    defaultValue="GSH",
)

stringGroup = VariableGroup(
    id="string_input",
    name="Input String",
    description="The input are in string",
    variables=[conservedResidues, glideOutputVariable, resNameProt, resNameLig, ligandName],
)
atomGroup = VariableGroup(
    id="atom_input",
    name="Input Atom",
    description="The input are in atom",
    variables=[conservedResidues, glideOutputVariable, residueProtein, residueLigand],
)

# Output variables
outputModelsVariable = PluginVariable(
    id="models",
    name="Alphafold models",
    description="The output models",
    type=VariableTypes.FOLDER,
)


# ==========================#
# Variable
# ==========================#
metricsVar = PluginVariable(
    name="Metrics ",
    id="metrics",
    description="The metrics list",
    type=VariableTypes.STRING,
    defaultValue="SG_S",
)


def finalAction(block: PluginBlock):

    bsc_result = block.inputs.get(glideOutputVariable.id, None)
    folder_to_analyse = bsc_result["dock_folder"]
    model_folder = bsc_result["model_folder"]

    conserved_indexes = block.inputs.get(conservedResidues.id, None)

    metrics = block.variables.get("metrics", "SG_S")

    if block.selectedInputGroup == stringGroup.id:
        res_name_prot = block.inputs.get(resNameProt.id, "CYS")
        res_name_lig = block.inputs.get(resNameLig.id, "SG")
        ligand_name = block.inputs.get(ligandName.id, "GSH")
    else:
        residue_protein = block.inputs.get(residueProtein.id, None)
        residue_ligand = block.inputs.get(residueLigand.id, None)

    models = prepare_proteins.proteinModels(model_folder)

    center_atom = {}  # Create dictionary to store the atom 3-element tuple for each model
    for model in models:  # Iterate the models inside the library
        # Iterate the residues for each Bio.PDB.Structure object
        for r in models.structures[model].get_residues():
            # Check that the residue matches the defined index
            aa = conserved_indexes[model]
            # for cons_ind in conserved_indexes[model]:
            if r.id[1] in conserved_indexes[model]:
                # Assert that the residue has the correct residue identity
                if r.resname == res_name_prot:
                    # Store the corresponsing tuple.
                    center_atom[model] = (r.get_parent().id, r.id[1], res_name_lig)
                    break

    print(f"center_atom: {center_atom}")

    atom_pairs = {}  # Define the dictionary containing the atom pairs for each model
    for model in models:
        atom_pairs[model] = {}
        for ligand in [ligand_name]:
            atom_pairs[model][ligand] = []
            atom_pairs[model][ligand].append((center_atom[model], "S1"))

    print(f"Atom pairs: {atom_pairs}")

    models.analyseDocking(folder_to_analyse, atom_pairs=atom_pairs)

    metric_distances = {}  # Define the global dictionary
    metric_distances[metrics] = {}  # Define the metric nested dictionary
    for model in models:
        metric_distances[metrics][model] = {}  # Define the model nested dictionary
        for ligand in models.docking_ligands[model]:
            # Define the ligand nested dictionary with all the docking distances list
            metric_distances[metrics][model][ligand] = models.getDockingDistances(model, ligand)

    print(f"metric_distances: {metric_distances}")

    models.combineDockingDistancesIntoMetrics(metric_distances)

    print(f"models.docking_data: {models.docking_data}")

    best_poses = models.getBestDockingPosesIteratively(metric_distances)

    models.extractDockingPoses(best_poses, folder_to_analyse, "best_docking_poses", separator="@")

    block.setOutput(outputModelsVariable.id, "best_docking_poses")


AnalyseGPXBlock = PluginBlock(
    name="Analyse Glide GPX",
    description="To analyse Glide GPX results",
    action=finalAction,
    variables=[metricsVar],
    inputGroups=[stringGroup, atomGroup],
    outputs=[outputModelsVariable],
)
