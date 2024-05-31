from HorusAPI import PluginBlock, PluginVariable, VariableGroup, VariableTypes

# TODO  Configure the inputs correctly

# ==========================#
# Variable inputs
# ==========================#
glideOutputVariable = PluginVariable(
    id="glide_output",
    name="Glide output",
    description="Glide output from the BSC calculations block",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_results"],
)
conservedResidues = PluginVariable(
    name="Conserved residues",
    id="conserved_indexes",
    description="The conserved residues",
    type=VariableTypes.CUSTOM,
    defaultValue=None,
)
residueProtein = PluginVariable(
    name="Atom Protein",
    id="resi_id1",
    description="Atom of the protein to calculate the distance to",
    type=VariableTypes.ATOM,
)
residueLigand = PluginVariable(
    name="Atom Ligand",
    id="resi_id2",
    description="Atom of the ligand to calculate the distance to",
    type=VariableTypes.ATOM,
)
resNameProt = PluginVariable(
    name="Protein residue name",
    id="res_name_prot",
    description="The protein residue name",
    type=VariableTypes.STRING,
    defaultValue="CYS",
)
atomNameProt = PluginVariable(
    name="Protein atomname",
    id="atom_name_prot",
    description="The protein atom name",
    type=VariableTypes.STRING,
    defaultValue="SG",
)
ligandName = PluginVariable(
    name="Ligand name",
    id="ligand_name",
    description="The ligand name",
    type=VariableTypes.STRING,
    defaultValue="GSH",
)
atomNameLig = PluginVariable(
    name="Ligand atom name",
    id="atom_name_ligand",
    description="The atom name of the ligand",
    type=VariableTypes.STRING,
    defaultValue="S1",
)

stringGroup = VariableGroup(
    id="string_input",
    name="Input String",
    description="The input are in string",
    variables=[
        conservedResidues,
        glideOutputVariable,
        resNameProt,
        atomNameProt,
        ligandName,
        atomNameLig,
    ],
)
atomGroup = VariableGroup(
    id="atom_input",
    name="Input Atom",
    description="The input are in atom",
    variables=[conservedResidues, glideOutputVariable, residueProtein, residueLigand],
)

# Output variables
outputModelsVariable = PluginVariable(
    id="best_poses",
    name="Best poses",
    description="The best poses from the analysis",
    type=VariableTypes.FOLDER,
)
analyseGlideOutputVariable = PluginVariable(
    id="glide_results_output",
    name="Glide results output",
    description="Output results of the Glide analysis",
    type=VariableTypes.CUSTOM,
    allowedValues=["glide_output"],
)

# ==========================#
# Variable
# ==========================#
metricsVar = PluginVariable(
    name="Metrics ",
    id="metrics",
    description="The metrics list",
    type=VariableTypes.STRING,
    defaultValue="SG_S",
)
removePreviousVar = PluginVariable(
    name="Remove previous models",
    id="remove_previous",
    description="Remove previous",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
separatorVar = PluginVariable(
    name="Separator",
    id="separator",
    description="The separator",
    type=VariableTypes.STRING,
    defaultValue="@",
)


def finalAction(block: PluginBlock):

    import prepare_proteins

    bsc_result = block.inputs.get(glideOutputVariable.id, None)
    folder_to_analyse = bsc_result["dock_folder"]
    model_folder = bsc_result["model_folder"]

    conserved_indexes = block.inputs.get(conservedResidues.id, None)

    metrics = block.variables.get("metrics", "SG_S")
    remove_previous = block.variables.get("remove_previous", False)
    separator = block.variables.get("separator", "@")

    if block.selectedInputGroup == stringGroup.id:
        res_name_prot = block.inputs.get(resNameProt.id, "CYS")
        atom_name_prot = block.inputs.get(atomNameProt.id, "SG")
        ligand_name = block.inputs.get(ligandName.id, "GSH")
        atom_name_lig = block.inputs.get(atomNameLig.id, "S1")
    else:
        residue_protein = block.inputs.get(residueProtein.id, None)
        res_name_prot = residue_protein["auth_comp_id"]
        atom_name_prot = residue_protein["auth_atom_id"]
        residue_ligand = block.inputs.get(residueLigand.id, None)
        ligand_name = residue_ligand["auth_comp_id"]
        atom_name_lig = residue_ligand["auth_atom_id"]
        metrics = f"{atom_name_prot}_{atom_name_lig}"

    models = prepare_proteins.proteinModels(model_folder)

    if conserved_indexes is None:
        raise ValueError("Conserved residues must be provided")
    if not isinstance(conserved_indexes, dict):
        try:
            conserved_indexes = int(conserved_indexes)
        except ValueError:
            raise ValueError("Conserved indexes must be an integer or a dictionary of integers")
        conserved_indexes_f = {}
        for model in models:
            conserved_indexes_f[model] = [conserved_indexes]
        conserved_indexes = conserved_indexes_f

    center_atom = {}  # Create dictionary to store the atom 3-element tuple for each model
    for model in models:  # Iterate the models inside the library
        # Iterate the residues for each Bio.PDB.Structure object
        for r in models.structures[model].get_residues():
            # Check that the residue matches the defined index
            # for cons_ind in conserved_indexes[model]:
            if r.id[1] in conserved_indexes[model]:
                # Assert that the residue has the correct residue identity
                if r.resname == res_name_prot:
                    # Store the corresponsing tuple.
                    center_atom[model] = (r.get_parent().id, r.id[1], atom_name_prot)
                    break

    atom_pairs = {}  # Define the dictionary containing the atom pairs for each model
    for model in models:
        atom_pairs[model] = {}
        for ligand in [ligand_name]:
            atom_pairs[model][ligand] = []
            atom_pairs[model][ligand].append((center_atom[model], atom_name_lig))

    models.analyseDocking(folder_to_analyse, atom_pairs=atom_pairs, separator=separator)

    metric_distances = {}  # Define the global dictionary
    metric_distances[metrics] = {}  # Define the metric nested dictionary
    for model in models:
        metric_distances[metrics][model] = {}  # Define the model nested dictionary
        for ligand in models.docking_ligands[model]:
            # Define the ligand nested dictionary with all the docking distances list
            metric_distances[metrics][model][ligand] = models.getDockingDistances(model, ligand)

    models.combineDockingDistancesIntoMetrics(metric_distances)

    best_poses = models.getBestDockingPosesIteratively(metric_distances)

    models.extractDockingPoses(
        best_poses,
        folder_to_analyse,
        "best_docking_poses",
        separator=separator,
        remove_previous=remove_previous,
    )

    block.setOutput(outputModelsVariable.id, "best_docking_poses")

    glideOutput = {
        "poses_folder": "best_docking_poses",
        "models_folder": model_folder,  # "prepared_proteins",
        "atom_pairs": atom_pairs,
    }
    import pickle

    with open("glide_output.pkl", "wb") as f:
        pickle.dump(glideOutput, f)

    block.setOutput(analyseGlideOutputVariable.id, glideOutput)


AnalyseGBlock = PluginBlock(
    name="Analyse Glide",
    description="To analyse Glide results",
    action=finalAction,
    variables=[metricsVar, removePreviousVar, separatorVar],
    inputGroups=[atomGroup, stringGroup],
    outputs=[outputModelsVariable, analyseGlideOutputVariable],
)
