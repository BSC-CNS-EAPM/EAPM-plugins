"""
Module containing the Mafft block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
proteinFolderVariable = PluginVariable(
    id="protein_folder",
    name="Protein folder",
    description="Folder containing the proteins",
    type=VariableTypes.FOLDER,
    defaultValue="proteins",
)

# ==========================#
# Variable outputs
# ==========================#
msaFile = PluginVariable(
    id="msa_file",
    name="MSA File",
    description="Multiple Sequence Alignment file",
    type=VariableTypes.FILE,
    allowedValues=["fasta"],
)


def calculate_msa(block: PluginBlock):
    """
    Calculate the Multiple Sequence Alignment (MSA) using Mafft.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """

    # pylint: disable=import-outside-toplevel
    import os
    import subprocess

    import bioprospecting
    import prepare_proteins

    # pylint: enable=import-outside-toplevel
    # Check that there is at least one pdb file in the folder
    protein_folder = block.inputs.get(proteinFolderVariable.id, "proteins")

    has_pdb = False
    for file in os.listdir(protein_folder):
        if file.endswith(".pdb"):
            has_pdb = True
            break

    if not has_pdb:
        raise ValueError(f"There are no pdb files in the protein folder: {protein_folder}")

    models = prepare_proteins.proteinModels(protein_folder)

    old_subprocess = subprocess.run(check=True)

    if block.remote.name != "Local":
        raise ValueError("This block only works on the local machine.")

    mafft_executable = block.config.get("mafft_path", "mafft")

    def hook_subprocess_mafft(command, **kwargs):
        if command.startswith("mafft"):
            command = command.replace("mafft", mafft_executable)
        return old_subprocess(command, **kwargs)

    try:
        subprocess.run = hook_subprocess_mafft
        msa = models.calculateMSA()

        output = "output.maf"

        bioprospecting.alignment.mafft.writeMSAToFile(msa, output, "clustal")

        block.setOutput(msaFile.id, output)
    except ValueError as e:
        raise ValueError(f"Error running Mafft: {e}") from e
    finally:
        subprocess.run = old_subprocess


multipleSequenceAlignmentBlock = PluginBlock(
    name="Multiple Sequence Alignment with Mafft",
    id="Mafft",
    description="Get the MSA from Mafft",
    inputs=[proteinFolderVariable],
    variables=[],
    outputs=[msaFile],
    action=calculate_msa,
)
