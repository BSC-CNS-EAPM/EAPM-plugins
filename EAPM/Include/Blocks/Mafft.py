"""
Module containing the Mafft block for the EAPM plugin
"""

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
msaVariable = PluginVariable(
    id="msa",
    name="MSA",
    description="Multiple Sequence Alignment",
    type=VariableTypes.CUSTOM,
    allowedValues=["msa"],
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
        block.setOutput("msaVariable", msa)
    finally:
        subprocess.run = oldSubprocess


multipleSequenceAlignmentBlock = PluginBlock(
    name="MultipleSequenceAlignment with Mafft",
    description="Get the MSA from Mafft",
    inputs=[proteinFolderVariable],
    variables=[],
    outputs=[msaVariable],
    action=calculateMSAAction,
)



