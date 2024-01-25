from HorusAPI import PluginBlock, PluginVariable, VariableTypes, Extensions

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


def getConservedMSAPositions(block: PluginBlock):
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
        conserved_index = models.getConservedMSAPositions(msa)

        # Generate an HTML with a table showing the conserved residues
        import pandas as pd

        df = pd.DataFrame(conserved_index, columns=["Index", "Residue"])

        html = df.to_html(index=False)

        Extensions().loadHTML(html, title="Conserved residues")

        # Get the residue index to get
        residueIndexes = block.variables.get("residue_index", [])

        if residueIndexes is None or len(residueIndexes) == 0:
            # Get all the indexes
            residueIndexes = []
            for index, _ in conserved_index:
                residueIndexes.append(index)
        else:
            # Convert the indexes to int
            residueIndexes = [int(index) for index in residueIndexes]

        conservedResidues = models.getStructurePositionsFromMSAindexes(residueIndexes)

        for model in conservedResidues:
            if len(conservedResidues[model]) == 0:
                raise Exception(
                    "There are no conserved residues for the selected indexes: "
                    + " ".join(residueIndexes)
                )
            break

        print(f"Conserved residues with index {residueIndexes}:", conservedResidues)

        block.setOutput("common_residues", conservedResidues)
    finally:
        subprocess.run = oldSubprocess


multipleSequenceAlignmentBlock = PluginBlock(
    name="Multiple sequence alignment",
    description="Get the conserved residues from a set of proteins",
    inputs=[proteinFolderVariable],
    variables=[residueIndexToGetVariable],
    outputs=[commonResiduesOutputVariable],
    action=getConservedMSAPositions,
)
