"""
Module containing the Align block for the EAPM plugin (mafft needed)
"""

from HorusAPI import PluginVariable, VariableTypes, PluginBlock, VariableList

# ==========================#
# Variable inputs
# ==========================#
inputFolderAlign = PluginVariable(
    name="Input folder",
    id="input_folder",
    description="The input folder with the PDBs to align.",
    type=VariableTypes.FOLDER,
    defaultValue="trimmed_models",
)
pdbReferenceAlign = PluginVariable(
    name="PDB reference",
    id="pdb_reference",
    description="The reference PDB to align to.",
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
integerChainIndexVariable = PluginVariable(
    name="Chain index",
    id="chain_index",
    description="Chain index.",
    type=VariableTypes.INTEGER,
    defaultValue=0,
)

chainIndexesAlign = VariableList(
    name="Chain indexes",
    id="chain_indexes",
    description="Chain indexes to use for the alignment. Use this option when the trajectories have corresponding chains in their topologies.",
    prototypes=[integerChainIndexVariable],
)

##############################
# Block's advanced variables #
##############################
trajectoryChainIndexVariable = PluginVariable(
    name="Trajectory chain index",
    id="trajectory_chain_index",
    description="Chain index number.",
    type=VariableTypes.INTEGER,
    defaultValue=0,
)

trajectoryChainIndexesAlign = VariableList(
    name="Trajectory chain indexes",
    id="trajectory_chain_indexes",
    description="Chain indexes of the target trajectories to use in the alignment.",
    prototypes=[trajectoryChainIndexVariable],
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
    name="Reference residue index",
    id="reference_residues",
    description="Reference residues.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
referenceResiduesAlign = VariableList(
    name="Reference residues",
    id="reference_residues",
    description="Reference residues.",
    prototypes=[referenceResiduesAlign],
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
    outputFolder = block.outputs.get("path", "aligned_models")
    chainIndexes = block.variables.get("chain_indexes", None)
    trajectoryChainIndexes = block.variables.get("trajectory_chain_indexes", None)
    alignmentMode = block.variables.get("alignment_mode", None)
    referenceResidues = block.variables.get("reference_residues", None)

    # Parse the chain indexes
    chainIndexes = [x["chain_index"] for x in chainIndexes]

    # Parse the trajectory chain indexes
    trajectoryChainIndexes = [x["trajectory_chain_index"] for x in trajectoryChainIndexes]

    # Parse the reference residues
    referenceResidues = [x["reference_residues"] for x in referenceResidues]

    import prepare_proteins

    print("Loading PDB files...")

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
    inputs=[pdbReferenceAlign, inputFolderAlign],
    outputs=[outputAlign],
)
