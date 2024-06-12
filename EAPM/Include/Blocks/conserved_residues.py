"""
Module containing the conserved residues block for the EAPM plugin
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

proteinFolderVariable = PluginVariable(
    id="protein_folder",
    name="Protein folder",
    description="Folder containing the proteins",
    type=VariableTypes.FOLDER,
    defaultValue="proteins",
)

residueIndexToGetVariable = PluginVariable(
    id="residue_index",
    name="Residue index",
    description="Index of the residues to get. If not set all the residues will be returned.",
    type=VariableTypes.LIST,
)

commonResiduesOutputVariable = PluginVariable(
    id="common_residues",
    name="Common residues",
    description="Common residues between the proteins",
    type=VariableTypes.CUSTOM,
    allowedValues=["multimodel_common_residue_grid"],
)


def get_conserved(block: PluginBlock):
    """
    Get the conserved residues from a set of proteins.

    Args:
        block (PluginBlock): The PluginBlock object.

    Returns:
        None
    """
    # pylint: disable=import-outside-toplevel
    import os
    import subprocess

    import pandas as pd
    import prepare_proteins

    # pylint: enable=import-outside-toplevel

    protein_folder = block.inputs.get(proteinFolderVariable.id, "proteins")

    # Check that there is at least one pdb file in the folder

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
        conserved_index = models.getConservedMSAPositions(msa)

        # Generate an HTML with a table showing the conserved residues

        df = pd.DataFrame(conserved_index, columns=["Index", "Residue"])

        html = df.to_html(index=False)

        Extensions().loadHTML(html, title="Conserved residues")

        # Get the residue index to get
        residue_indexes = block.variables.get(residueIndexToGetVariable.id, [])

        if residue_indexes is None or len(residue_indexes) == 0:
            # Get all the indexes
            residue_indexes = []
            for index, _ in conserved_index:
                residue_indexes.append(index)
        else:
            # Convert the indexes to int
            residue_indexes = [int(index) for index in residue_indexes]

        conserved_residues = models.getStructurePositionsFromMSAindexes(residue_indexes)

        for model in conserved_residues:
            if len(conserved_residues[model]) == 0:
                raise ValueError(
                    "There are no conserved residues for the selected indexes: "
                    + " ".join(str(residue_indexes))
                )
            break

        print(f"Conserved residues with index {residue_indexes}:", conserved_residues)

        block.setOutput("common_residues", conserved_residues)
    finally:
        subprocess.run = old_subprocess


conservedResiduesMSABlock = PluginBlock(
    name="Conserved Residues from MSA",
    id="conserved_residues_msa",
    description="Get the conserved residues from a set of proteins",
    inputs=[proteinFolderVariable],
    variables=[residueIndexToGetVariable],
    outputs=[commonResiduesOutputVariable],
    action=get_conserved,
)
