"""
Module containing the Mafft block for the EAPM plugin
"""
import Bio.AlignIO, Bio.SeqIO
from HorusAPI import PluginBlock, PluginVariable, VariableTypes, Extensions


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


def calculateMSAAction(block: PluginBlock):
    proteinFolder = block.inputs.get("protein_folder", "proteins")

    import prepare_proteins

    # Check that there is at least one pdb file in the folder
    import os

    hasPDB = False
    for file in os.listdir(proteinFolder):
        if file.endswith(".pdb"):
            hasPDB = True
            break

    if not hasPDB:
        raise Exception(f"There are no pdb files in the protein folder: {proteinFolder}")

    models = prepare_proteins.proteinModels(proteinFolder)

    import subprocess

    oldSubprocess = subprocess.run

    if block.remote.name != "Local":
        raise Exception("This block only works on the local machine.")

    mafftExecutable = block.config.get("mffa_path", "mffa")

    def hookSubprocessMafft(command, **kwargs):
        if command.startswith("mafft"):
            command = command.replace("mafft", mafftExecutable)
        return oldSubprocess(command, **kwargs)

    try:
        subprocess.run = hookSubprocessMafft
        msa = models.calculateMSA()
        
        output = block.outputs.get("msa_file", "output.fasta")
        with open(output, "w") as output_handle:
            AlignIO.write(msa, output_handle, "fasta")
            
        block.setOutput("msaVariable", output)
    finally:
        subprocess.run = oldSubprocess


multipleSequenceAlignmentBlock = PluginBlock(
    name="Multiple Sequence Alignment with Mafft",
    description="Get the MSA from Mafft",
    inputs=[proteinFolderVariable],
    variables=[],
    outputs=[msaFile],
    action=calculateMSAAction,
)



