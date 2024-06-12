"""
Module containing the analyse PELE block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

peleOutputFolderInput = PluginVariable(
    id="pele_folder",
    name="Pele folder",
    description="Folder containing PELE output",
    type=VariableTypes.FOLDER,
    defaultValue="pele",
)


def analyse_PELE(block: PluginBlock):
    """
    Analyze PELE data and calculates catalytic distances.

    Args:
        block (PluginBlock): The PluginBlock object representing the PELE block.

    Returns:
        None

    Raises:
        None
    """
    # pylint: disable=import-outside-toplevel
    import pele_analysis

    # pylint: enable=import-outside-toplevel

    pele_folder = block.inputs.get(peleOutputFolderInput.id, "pele")

    pele = pele_analysis.peleAnalysis(
        pele_folder,
        verbose=True,
        separator="-",
        trajectories=False,
        data_folder_name="pele_data/",
        read_equilibration=True,
    )

    # Classify distances into common metrics --> the catalytic labels
    catalytic_names = ["SER-L", "SER-HIS", "HIS-ASP"]

    catalytic_labels = {}
    for cn in catalytic_names:
        catalytic_labels[cn] = {}
        for protein in pele.proteins:
            if protein not in catalytic_labels[cn]:
                catalytic_labels[cn][protein] = {}
            for ligand in pele.ligands:
                if ligand not in catalytic_labels[cn][protein]:
                    catalytic_labels[cn][protein][ligand] = []

    for protein in pele.proteins:
        for ligand in pele.ligands:
            distances = pele.getDistances(protein, ligand)
            if distances is None:
                continue
            for d in distances:
                at1 = d.replace("distance_", "").split("_")[0]
                at2 = d.replace("distance_", "").split("_")[1]
                if at1.endswith("OG") and at2.startswith("L"):
                    catalytic_labels["SER-L"][protein][ligand].append(d)

                elif at1.endswith("OG") and at2.endswith("NE2"):
                    catalytic_labels["SER-HIS"][protein][ligand].append(d)

                elif at1.endswith("ND1") and (at2.endswith("OD1") or at2.endswith("OD2")):
                    catalytic_labels["HIS-ASP"][protein][ligand].append(d)

    print("catalytic_labels", catalytic_labels)

    # Calculate the catalytic distances for each group of distances
    pele.combineDistancesIntoMetrics(catalytic_labels, overwrite=True)

    print("Success combining distances into metrics")

    pele.bindingEnergyLandscape()


analysePELEBlock = PluginBlock(
    name="Analyse PELE",
    id="analyse_pele",
    description="Analyse PELE output",
    action=analyse_PELE,
    inputs=[peleOutputFolderInput],
    variables=[],
    outputs=[],
)
