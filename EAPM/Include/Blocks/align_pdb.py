"""
Module containing the Align block for the EAPM plugin (mafft needed)
"""

from HorusAPI import PluginBlock, PluginVariable, VariableList, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
inputFolder = PluginVariable(
    name="Input folder",
    id="input_folder",
    description="The input folder with the PDBs to align.",
    type=VariableTypes.FOLDER,
    defaultValue="trimmed_models",
)
pdbReference = PluginVariable(
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
output = PluginVariable(
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
chainIndexes = VariableList(
    name="Chain indexes",
    id="chain_indexes",
    description="Chain indexes to use for the alignment. "
    "Use this option when the trajectories have corresponding chains in their topologies.",
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
trajectoryChainIndexes = VariableList(
    name="Trajectory chain indexes",
    id="trajectory_chain_indexes",
    description="Chain indexes of the target trajectories to use in the alignment.",
    prototypes=[trajectoryChainIndexVariable],
)
alignmentMode = PluginVariable(
    name="Alignment mode",
    id="alignment_mode",
    description="The mode defines how sequences are aligned. "
    "'exact' for structurally aligning positions with exactly the same amino acids after the sequence alignment "
    "or 'aligned' for structurally aligning sequences using all positions aligned in the sequence alignment.",
    type=VariableTypes.STRING_LIST,
    defaultValue="aligned",
    allowedValues=["aligned", "exact"],
)
integerReferenceResidues = PluginVariable(
    name="Reference residue index",
    id="reference_residues_int",
    description="Reference residues.",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
referenceResidues = VariableList(
    name="Reference residues",
    id="reference_residues",
    description="Reference residues.",
    prototypes=[integerReferenceResidues],
)


# Align action block
def initial_action(block: PluginBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    # pylint: disable=import-outside-toplevel
    import subprocess

    import prepare_proteins

    # pylint: enable=import-outside-toplevel

    if block.remote.name != "Local":
        raise ValueError("This block is only available for local execution.")

    mafft_executable = block.config.get("mafft_path", "mafft")
    # Check if the remote has the mafft executable
    try:
        block.remote.remoteCommand(f"{mafft_executable} --version")
    except Exception as exc:
        raise ValueError(
            "MAFFT is not installed in the selected machine. "
            "Please install MAFFT in order to align the PDBs"
        ) from exc

    # Loading plugin variables
    input_folder = block.inputs.get(inputFolder.id, None)
    pdb_reference = block.inputs.get(pdbReference.id, None)
    output_folder = block.outputs.get(output.id, "aligned_models")

    default_chain_index = {"chain_index": 0}

    chain_indexes = block.variables.get(chainIndexes.id, [default_chain_index])
    trajectory_chain_indexes = block.variables.get(trajectoryChainIndexes.id, [])
    alignment_mode = block.variables.get(alignmentMode.id, "aligned")
    reference_residues = block.variables.get(referenceResidues.id, [])

    print("Loading PDB files...")

    models = prepare_proteins.proteinModels(input_folder)

    # Parse the chain indexes
    if chain_indexes is not None:
        chain_indexes = [x["chain_index"] for x in chain_indexes]
    else:
        chain_indexes = [0]

    trajectory_chain_indexes = None
    # Parse the trajectory chain indexes
    if trajectory_chain_indexes is not None:
        trajectory_chain_indexes = [x["trajectory_chain_index"] for x in trajectory_chain_indexes]
        trajectory_chain_indexes = {}
        for i, model in enumerate(models.models_names):
            trajectory_chain_indexes[model] = trajectory_chain_indexes[i]

    # Parse the reference residues
    if reference_residues is not None:
        reference_residues = [x["reference_residues"] for x in reference_residues]

    print("Aligning models...")

    old_subprocess = subprocess.run

    def hook_subprocess_mafft(command, **kwargs):
        if command.startswith("mafft"):
            command = command.replace("mafft", mafft_executable)

        print("Running command:", command)
        return old_subprocess(command, check=True, **kwargs)

    try:
        subprocess.run = hook_subprocess_mafft
        models.alignModelsToReferencePDB(
            pdb_reference,
            output_folder,
            chain_indexes=chain_indexes,
            trajectory_chain_indexes=trajectory_chain_indexes,
            aligment_mode=alignment_mode,
            reference_residues=reference_residues,
            verbose=True,
        )
    finally:
        subprocess.run = old_subprocess

    print("Setting output of block to the results directory...")

    # Set the output
    block.setOutput(output.id, output_folder)


alignBlock = PluginBlock(
    name="Align PDBs",
    id="align_PDBs",
    description="Align all models to a reference PDB based on a sequence alignment. (For local)",
    action=initial_action,
    variables=[
        chainIndexes,
        trajectoryChainIndexes,
        alignmentMode,
        referenceResidues,
    ],
    inputs=[pdbReference, inputFolder],
    outputs=[output],
)
