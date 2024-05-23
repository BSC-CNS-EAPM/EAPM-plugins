import json
import os
import shutil

import bsc_calculations
import pandas as pd
import prepare_proteins

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
fasta_fileAF = PluginVariable(
    name="Fasta file",
    id="fasta_file",
    description="The input fasta file.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# Output variables
outputModelsVariable = PluginVariable(
    id="models",
    name="Alphafold models",
    description="The output models",
    type=VariableTypes.FOLDER,
)


def finalAction(block: PluginBlock):

    models = prepare_proteins.proteinModels("models")

    with open("lig_atom_name.json", "r") as f:
        lig_atom_name = json.load(f)

    triads = {}
    triads["FeLip9"] = [111, 190, 167]

    triad_atoms = {}
    triad_atoms["ser_OG"] = {}
    triad_atoms["his_ND1"] = {}
    triad_atoms["his_NE2"] = {}
    triad_atoms["asp_OD1"] = {}
    triad_atoms["asp_OD2"] = {}

    for model in models:  # Iterate the models inside the library
        S = triads[model][0]
        H = triads[model][1]
        D = triads[model][2]
        for r in models.structures[
            model
        ].get_residues():  # Iterate the residues for each Bio.PDB.Structure object
            if r.id[1] == S:  # Check that the residue matches the defined index
                assert (
                    r.resname == "SER"
                )  # Assert that the residue has the correct residue identity
                triad_atoms["ser_OG"][model] = (
                    r.get_parent().id,
                    r.id[1],
                    "OG",
                )  # Store the corresponsing tuple.
            elif r.id[1] == H:
                assert r.resname == "HIS"
                triad_atoms["his_ND1"][model] = (r.get_parent().id, r.id[1], "ND1")
                triad_atoms["his_NE2"][model] = (r.get_parent().id, r.id[1], "NE2")
            elif r.id[1] == D:
                assert r.resname == "ASP"
                triad_atoms["asp_OD1"][model] = (r.get_parent().id, r.id[1], "OD1")
                triad_atoms["asp_OD2"][model] = (r.get_parent().id, r.id[1], "OD2")

    atom_pairs = {}  # Define the dictionary containing the atom pairs for each model
    for model in models:
        atom_pairs[model] = {}
        atom_pairs[model]["PET"] = []
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[2]["C1"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[2]["C8"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[4]["C1"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[4]["C8"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[6]["C1"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[6]["C8"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[8]["C1"]))
        atom_pairs[model]["PET"].append((triad_atoms["ser_OG"][model], lig_atom_name[8]["C8"]))

    models.analyseDocking("docking", atom_pairs=atom_pairs)

    metric_distances = {}  # Define the global dictionary
    metric_distances["OG_C"] = {}  # Define the metric nested dictionary

    for model in models:
        metric_distances["OG_C"][model] = {}  # Define the model nested dictionary
        for ligand in models.docking_ligands[model]:
            # Define the ligand nested dictionary with all the docking distances list
            metric_distances["OG_C"][model][ligand] = models.getDockingDistances(model, ligand)

    models.combineDockingDistancesIntoMetrics(metric_distances)

    best_poses = models.getBestDockingPosesIteratively(metric_distances)

    models.extractDockingPoses(best_poses, "docking", "best_docking_poses")

    block.setOutput(outputModelsVariable.id, "best_docking_poses")
