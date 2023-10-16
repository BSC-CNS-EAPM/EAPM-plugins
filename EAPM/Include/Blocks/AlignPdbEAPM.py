"""
Module containing the AlphaFold block for the NBDSuite plugin
"""

import datetime
import os
import shutil
import subprocess
import tarfile
from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

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
    name="Alphafold output",
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
    description=" Chain indexes to use for the alignment. Use this option when the trajectories have corresponding chains in their topologies.",
    type=VariableTypes.ANY,
    defaultValue=None,
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
    defaultValue="align",
    allowedValues=["align", "exact"],
)
referenceResiduesAlign = PluginVariable(
    name="Reference residues",
    id="reference_residues",
    description="Reference residues.",
    type=VariableTypes.STRING_LIST,
    defaultValue=None,
)


# Align action block
def initialAlign(block: SlurmBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    # Loading plugin variables
    inputFolder = block.inputs.get("input_folder", "None")
    pdbReference = block.inputs.get("pdb_reference", "None")
    outputFolder = block.inputs.get("path", "None")
    chainIndexes = block.inputs.get("chain_indexes", "None")
    trajectoryChainIndexes = block.variables.get("trajectory_chain_indexes", "None")
    alignmentMode = block.variables.get("alignment_mode", "None")
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
        alignment_mode=alignmentMode,
        reference_residues=referenceResidues,
    )


# Block's final action
def finalAlign(block: SlurmBlock):
    """
    Final action of the block. It downloads the results from the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    outputFolder = block.inputs.get("path", "None")
    print("Setting output of block to the results directory...")

    # Set the output
    block.setOutput("path", outputFolder)


alignBlock = SlurmBlock(
    name="Align PDBs",
    description="Align all models to a reference PDB based on a sequence alignment. (For local)",
    initialAction=initialAlign,
    finalAction=finalAlign,
    variables=[
        chainIndexesAlign,
        trajectoryChainIndexesAlign,
        alignmentModeAlign,
        referenceResiduesAlign,
    ],
    inputs=[inputFolderAlign, pdbReferenceAlign],
    outputs=[outputAlign],
)
