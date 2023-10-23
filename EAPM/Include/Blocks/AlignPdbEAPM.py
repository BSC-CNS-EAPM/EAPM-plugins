"""
Module containing the Align block for the EAPM plugin (mafft needed)
"""

import datetime
import os
import shutil
import subprocess
import tarfile
from HorusAPI import PluginVariable, VariableTypes, PluginBlock

# ==========================#
# Variable inputs
# ==========================#
inputFolderAlign = PluginVariable(
    name="Input folder",
    id="input_folder",
    description="The input folder with the Pdbs to align.",
    type=VariableTypes.FOLDER,
    defaultValue="trimmed_models",
)
pdbReferenceAlign = PluginVariable(
    name="Pdb reference",
    id="pdb_reference",
    description="The reference Pdb to align to.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["pdb"],
)

# ==========================#
# Variable outputs
# ==========================#
outputAlign = PluginVariable(
    name="Align output",
    id="path",
    description="The folder containing the results.",
    type=VariableTypes.FOLDER,
    defaultValue="aligned_models",
)

##############################
#       Other variables      #
##############################
chainIndexesAlign = PluginVariable(
    name="Chain indexes",
    id="chain_indexes",
    description="Chain indexes to use for the alignment. Use this option when the trajectories have corresponding chains in their topologies.",
    type=VariableTypes.INTEGER_LIST,
    defaultValue=0,
)

##############################
# Block's advanced variables #
##############################
trajectoryChainIndexesAlign = PluginVariable(
    name="Trajectory chain indexes",
    id="trajectory_chain_indexes",
    description="Chain indexes of the target trajectory to use in the alignment.",
    type=VariableTypes.INTEGER_LIST,
    defaultValue=None,
)
alignmentModeAlign = PluginVariable(
    name="Alignment mode",
    id="alignment_mode",
    description="The mode defines how sequences are aligned. 'exact' for structurally aligning positions with exactly the same aminoacids after the sequence alignment or 'aligned' for structurally aligning sequences using all positions aligned in the sequence alignment.",
    type=VariableTypes.STRING_LIST,
    defaultValue="aligned",
    allowedValues=["aligned", "exact"],
)
referenceResiduesAlign = PluginVariable(
    name="Reference residues",
    id="reference_residues",
    description="Reference residues.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
)


# Align action block
def initialAlign(block: PluginBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    # Loading plugin variables
    inputFolder = block.inputs.get("input_folder", "None")
    pdbReference = block.inputs.get("pdb_reference", "None")
    outputFolder = block.inputs.get("path", "aligned_models")
    chainIndexes = block.inputs.get("chain_indexes", 0)
    trajectoryChainIndexes = block.variables.get("trajectory_chain_indexes", "None")
    alignmentMode = block.variables.get("alignment_mode", "aligned")
    referenceResidues = block.variables.get("reference_residues", "None")

    import prepare_proteins

    print("Loading Pdbs files...")

    models = prepare_proteins.proteinModels(inputFolder)

    print("Aligning models...")

    models.alignModelsToReferencePDB(
        pdbReference,
        outputFolder,
        chain_indexes=chainIndexes,
        trajectory_chain_indexes=trajectoryChainIndexes,
        aligment_mode=alignmentMode,
        reference_residues=referenceResidues,
        verbose=True,
    )

    print("Setting output of block to the results directory...")

    # Set the output
    block.setOutput("path", outputFolder)


alignBlock = PluginBlock(
    name="Align PDBs",
    description="Align all models to a reference PDB based on a sequence alignment. (For local)",
    action=initialAlign,
    variables=[
        chainIndexesAlign,
        trajectoryChainIndexesAlign,
        alignmentModeAlign,
        referenceResiduesAlign,
    ],
    inputs=[inputFolderAlign, pdbReferenceAlign],
    outputs=[outputAlign],
)
