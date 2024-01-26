import random

from HorusAPI import SlurmBlock, PluginVariable, VariableTypes, VariableGroup, VariableList


# Input variables
yamlPELEFileVariable = PluginVariable(
    id="yaml_pele_file",
    name="PELE yaml",
    description="YAML file containing the PELE configuration",
    type=VariableTypes.FILE,
    defaultValue="cst_input.yaml",
    allowedValues=["yaml"],
)

modelFolderVariable = PluginVariable(
    id="model_folder",
    name="Model folder",
    description="Folder containing the models",
    type=VariableTypes.FOLDER,
    defaultValue="models",
)

bestDockingPosesVariable = PluginVariable(
    id="best_docking_poses",
    name="Best docking poses",
    description="Number of best docking poses to analyse",
    type=VariableTypes.INTEGER,
    defaultValue=100,
)

glideOutputVariable = PluginVariable(
    id="glide_output",
    name="Glide output",
    description="Glide output from the Glide docking block",
    type=VariableTypes.CUSTOM,
    allowedValues=["glide_output"],
)


folderInputGroup = VariableGroup(
    id="folder_input_group",
    name="Folder input group",
    description="Input the model and ligand folders after a Dcoking Grid setup has been run",
    variables=[modelFolderVariable, bestDockingPosesVariable, yamlPELEFileVariable],
)

glideOutputGroup = VariableGroup(
    id="glide_output_group",
    name="Glide output group",
    description="Input the Glide output from the Glide docking block",
    variables=[glideOutputVariable, yamlPELEFileVariable],
)

# Variables
peleFolderNameVariable = PluginVariable(
    id="pele_folder_name",
    name="PELE folder name",
    description="Name of the PELE folder",
    type=VariableTypes.STRING,
    defaultValue="pele",
)

atom1SelectionVariable = PluginVariable(
    id="protein_atom",
    name="Protein atom",
    description="Select an atom on the protein to calculate the distance",
    type=VariableTypes.ATOM,
)

atom2SelectionVariable = PluginVariable(
    id="ligand_atom",
    name="Ligand atom",
    description="Select an atom on the ligand to calculate the distance",
    type=VariableTypes.ATOM,
)

atomSelectionsGroupVariable = PluginVariable(
    id="group",
    name="Group",
    description="Name of the group selection to clusterize",
    type=VariableTypes.STRING,
)

overrideLigandNameVariable = PluginVariable(
    id="override_ligand_name",
    name="Override ligand name",
    description="Override the ligand name with a custom name",
    type=VariableTypes.STRING,
)

selectionsListVariable = VariableList(
    id="selections_list",
    name="Selections",
    description="List of selections to analyse",
    prototypes=[
        atom1SelectionVariable,
        atom2SelectionVariable,
        atomSelectionsGroupVariable,
        overrideLigandNameVariable,
    ],
)

# PELE variables
boxCentersVariable = PluginVariable(
    id="box_centers",
    name="Box centers",
    description="List of box centers",
    type=VariableTypes.SPHERE,
    category="PELE",
)

constraintsVariable = PluginVariable(
    id="constraints",
    name="Constraints",
    description="List of constraints (not implemented yet)",
    type=VariableTypes.STRING,
    category="PELE",
)

ligandIndexVariable = PluginVariable(
    id="ligand_index",
    name="Ligand index",
    description="Ligand index",
    type=VariableTypes.INTEGER,
    defaultValue=1,
    category="PELE",
)

peleStepsVariable = PluginVariable(
    id="pele_steps",
    name="PELE steps",
    description="Number of PELE steps",
    type=VariableTypes.INTEGER,
    defaultValue=100,
    category="PELE",
)

peleDebugVariable = PluginVariable(
    id="pele_debug",
    name="PELE debug",
    description="Debug mode",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

peleIterationsVariable = PluginVariable(
    id="pele_iterations",
    name="PELE iterations",
    description="Number of PELE iterations",
    type=VariableTypes.INTEGER,
    defaultValue=5,
    category="PELE",
)

equilibrationStepsVariable = PluginVariable(
    id="equilibration_steps",
    name="Equilibration steps",
    description="Number of equilibration steps",
    type=VariableTypes.INTEGER,
    defaultValue=100,
    category="PELE",
)

ligandEnergyGroupsVariable = PluginVariable(
    id="ligand_energy_groups",
    name="Ligand energy groups",
    description="List of ligand energy groups (not implemented yet)",
    type=VariableTypes.STRING,
    category="PELE",
)

peleSeparatorVariable = PluginVariable(
    id="pele_separator",
    name="PELE separator",
    description="Separator for the PELE models and ligands",
    type=VariableTypes.STRING,
    defaultValue="-",
    category="PELE",
)

usePeleffyVariable = PluginVariable(
    id="use_peleffy",
    name="Use PELEffy",
    description="Use PELEffy to generate the ligand parameters",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="PELE",
)

useSrunVariable = PluginVariable(
    id="use_srun",
    name="Use srun",
    description="Use srun to launch PELE",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="PELE",
)

energyByResidueVariable = PluginVariable(
    id="energy_by_residue",
    name="Energy by residue",
    description="Calculate the energy by residue",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

ebrNewFlagVariable = PluginVariable(
    id="ebr_new_flag",
    name="EBR new flag",
    description="New flag for the energy by residue calculation",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

ninetyDegreesVersionVariable = PluginVariable(
    id="ninety_degrees_version",
    name="90 degrees version",
    description="Use the 90 degrees version of PELE",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

analysisVariable = PluginVariable(
    id="analysis",
    name="Analysis",
    description="Perform analysis",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)


energyByResidueTypeVariable = PluginVariable(
    id="energy_by_residue_type",
    name="Energy by residue type",
    description="Type of energy by residue calculation",
    type=VariableTypes.STRING_LIST,
    allowedValues=["all", "lennard_jones", "sgb", "electrostatic"],
    defaultValue="all",
    category="PELE",
)

peptideVariable = PluginVariable(
    id="peptide",
    name="Peptide",
    description="Peptide",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

equilibrationModeVariable = PluginVariable(
    id="equilibration_mode",
    name="Equilibration mode",
    description="Equilibration mode",
    type=VariableTypes.STRING_LIST,
    allowedValues=["equilibrationLastSnapshot", "equilibrationCluster"],
    defaultValue="equilibrationLastSnapshot",
    category="PELE",
)

spawningVariable = PluginVariable(
    id="spawning",
    name="Spawning",
    description="Spawning",
    type=VariableTypes.STRING_LIST,
    allowedValues=[
        "independent",
        "inverselyProportional",
        "epsilon",
        "variableEpsilon",
        "independentMetric",
        "UCB",
        "FAST",
        "ProbabilityMSM",
        "MetastabilityMSM",
        "IndependentMSM",
    ],
    defaultValue="independent",
    category="PELE",
)

continuationVariable = PluginVariable(
    id="continuation",
    name="Continuation",
    description="Continuation",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

equilibrationVariable = PluginVariable(
    id="equilibration",
    name="Equilibration",
    description="Equilibration",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="PELE",
)

skipModelsVariable = PluginVariable(
    id="skip_models",
    name="Skip models",
    description="Write which model names to skip",
    type=VariableTypes.LIST,
    category="PELE",
)

skipLigandsVariable = PluginVariable(
    id="skip_ligands",
    name="Skip ligands",
    description="Write which ligands to skip",
    type=VariableTypes.LIST,
    category="PELE",
)

extendIterationsVariable = PluginVariable(
    id="extend_iterations",
    name="Extend iterations",
    description="Extend the number of iterations",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE",
)

logFileVariable = PluginVariable(
    id="log_file",
    name="Log file",
    description="Enable log file",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE"
)

rescoringVariable = PluginVariable(
    id="rescoring",
    name="Rescoring",
    description="Enable rescoring",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
    category="PELE"
)

epsilonVariable = PluginVariable(
    id="epsilon",
    name="Epsilon",
    description="TODO Epsilon description",
    type=VariableTypes.FLOAT,
    defaultValue=0.5,
    category="PELE"
)

ligandEquilibrationCstValue = PluginVariable(
    id="ligand_equilibration_cst",
    name="Ligand equilibration cst",
    description="TODO Ligand equilibration cst description",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="PELE"
)

onlyModelsVariable = PluginVariable(
    id="only_models",
    name="Only models",
    description="Write which model names to use",
    type=VariableTypes.LIST,
    category="PELE",
)

onlyLigandsVariable = PluginVariable(
    id="only_ligands",
    name="Only ligands",
    description="Write which ligands to use",
    type=VariableTypes.LIST,
    category="PELE",
)

onlyCombinationsVariable = PluginVariable(
    id="only_combinations",
    name="Only combinations",
    description="Write which combinations to use",
    type=VariableTypes.LIST,
    category="PELE",
)

ligandTemplateVariable = PluginVariable(
    id="ligand_template",
    name="Ligand template",
    description="Ligand template (not implemented yet)",
    type=VariableTypes.STRING,
    category="PELE",
)

seedVariable = PluginVariable(
    id="seed",
    name="Seed",
    description="Seed to use during the PELE calculation. If set to -1, a random seed will be used",
    type=VariableTypes.INTEGER,
    defaultValue=-1,
    category="PELE",
)

# Outputs
peleOutputFolderOutput = PluginVariable(
    id="pele_output_folder",
    name="PELE folder",
    description="Folder containing the PELE output",
    type=VariableTypes.FOLDER,
)


def peleAction(block: SlurmBlock):
    if block.selectedInputGroup == "glide_output_group":
        glide_output = block.inputs.get("glide_output")
        poses_folder = glide_output.get("poses_folder")
        models_folder = glide_output.get("models_folder")
        atom_pairs = glide_output.get("atom_pairs", {})

        if poses_folder is None or models_folder is None:
            raise Exception("Glide output does not contain the required folders")
    else:
        poses_folder = block.inputs.get("poses_folder")
        models_folder = block.inputs.get("models_folder")
        atom_pairs = {}

    # Get all the variables from the block
    boxCentersValue = block.variables.get("box_centers", [])
    constraintsValue = block.variables.get("constraints", [])
    ligandIndexValue = block.variables.get("ligand_index", 1)
    peleStepsValue = block.variables.get("pele_steps", 100)
    peleDebugValue = block.variables.get("pele_debug", False)
    peleIterationsValue = block.variables.get("pele_iterations", 5)
    equilibrationStepsValue = block.variables.get("equilibration_steps", 100)
    ligandEnergyGroupsValue = block.variables.get("ligand_energy_groups", [])
    peleSeparatorValue = block.variables.get("pele_separator", "-")
    usePeleffyValue = block.variables.get("use_peleffy", True)
    useSrunValue = block.variables.get("use_srun", True)
    energyByResidueValue = block.variables.get("energy_by_residue", False)
    ebrNewFlagValue = block.variables.get("ebr_new_flag", False)
    ninetyDegreesVersionValue = block.variables.get("ninety_degrees_version", False)
    analysisValue = block.variables.get("analysis", False)
    energyByResidueTypeValue = block.variables.get("energy_by_residue_type", "all")
    peptideValue = block.variables.get("peptide", False)
    equilibrationModeValue = block.variables.get(
        "equilibration_mode", "equilibrationLastSnapshot"
    )
    spawningValue = block.variables.get("spawning", "independent")
    continuationValue = block.variables.get("continuation", False)
    equilibrationValue = block.variables.get("equilibration", True)
    skipModelsValue = block.variables.get("skip_models", [])
    skipLigandsValue = block.variables.get("skip_ligands", [])
    extendIterationsValue = block.variables.get("extend_iterations", False)
    onlyModelsValue = block.variables.get("only_models", [])
    onlyLigandsValue = block.variables.get("only_ligands", [])
    onlyCombinationsValue = block.variables.get("only_combinations", [])
    ligandTemplateValue = block.variables.get("ligand_template", "")
    seedValue = block.variables.get("seed", -1)
    logFileValue = block.variables.get("log_file", False)
    rescoringValue = block.variables.get("rescoring", False)
    epsilonValue = block.variables.get("epsilon", 0.5)
    ligandEquilibrationCstValue = block.variables.get("ligand_equilibration_cst", True)

    # Parse spawningValue
    validSpawnings = ['independent', 'inverselyProportional', 'epsilon', 'variableEpsilon',
                     'independentMetric', 'UCB', 'FAST', 'ProbabilityMSM', 'MetastabilityMSM',
                     'IndependentMSM']
    
    if spawningValue != None and spawningValue not in validSpawnings:
            message = 'Spawning method %s not found.' % spawningValue
            message = 'Allowed options are: ' + str(validSpawnings)
            raise ValueError(message)

    # Parse energyByResidueValue
    energy_by_residue_types = ['all', 'lennard_jones', 'sgb', 'electrostatic']
    if energyByResidueValue not in energy_by_residue_types:
        raise ValueError('%s not found. Try: %s' % (energyByResidueValue, energy_by_residue_types))

    # Parse seedValue
    if seedValue == -1:
        seedValue = random.randint(0, 1000000)

    # Parse ligandEnergyGroups
    if not isinstance(ligandEnergyGroupsValue, type(None)):
        if not isinstance(ligandEnergyGroupsValue, dict):
            raise ValueError('Ligand energy groups, must be given as a dictionary')

    import prepare_proteins

    print("Using models folder: ", models_folder)
    models = prepare_proteins.proteinModels(models_folder)

    selections = block.variables.get("selections_list", [])
    if atom_pairs == {}:
        groups = []
        for model in models:
            atom_pairs[model] = {}
            for selection in selections:
                current_group = selection["group"]
                if current_group not in groups:
                    groups.append(current_group)
                atom1 = selection["protein_atom"]
                protein_chain = atom1["chainID"]
                prtoein_resnum = atom1["residue"]
                protein_atom = atom1["auth_atom_id"]
                protein_tuple = (protein_chain, prtoein_resnum, protein_atom)
                atom2 = selection["ligand_atom"]
                ligandName = atom2["auth_comp_id"]
                if (
                    selection.get("override_ligand_name") is not None
                    and selection["override_ligand_name"] != ""
                ):
                    ligandName = selection["override_ligand_name"]
                ligand_atom = atom2["auth_atom_id"]
                if atom_pairs[model].get(ligandName, None) is None:
                    atom_pairs[model][ligandName] = []
                atom_pairs[model][ligandName].append((protein_tuple, ligand_atom))

    cst_yaml = block.inputs.get("yaml_pele_file")

    cpus = block.variables.get("cpus", 48)
    peleFolderName = block.variables.get("pele_folder_name", "pele")

    # Remaining variables to implement:

    # box_centers=None, constraints=None, box_radius=10, ligand_energy_groups=None,
    # skip_models=None, skip_ligands=None,
    # only_models=None, only_ligands=None, only_combinations=None, ligand_templates=None,
    # nonbonded_energy=None, nonbonded_energy_type='all', nonbonded_new_flag=False, covalent_setup=False, covalent_base_aa=None,
    # membrane_residues=None, bias_to_point=None, com_bias1=None, com_bias2=None

    # Setup pele
    jobs = models.setUpPELECalculation(
        peleFolderName,
        poses_folder,
        cst_yaml,
        iterations=peleIterationsValue,
        cpus=cpus,
        distances=atom_pairs,
        separator=peleSeparatorValue,
        steps=peleStepsValue,
        seed=seedValue,
        energy_by_residue=energyByResidueValue,
        debug=peleDebugValue,
        equilibration_steps=equilibrationStepsValue,
        ligand_index=ligandIndexValue,
        use_peleffy=usePeleffyValue,
        usesrun=useSrunValue,
        ebr_new_flag=ebrNewFlagValue,
        ninety_degrees_version=ninetyDegreesVersionValue,
        extend_iterations=extendIterationsValue,
        continuation=continuationValue,
        equilibration=equilibrationValue,
        peptide=peptideValue,
        analysis=analysisValue,
        energy_by_residue_type=energyByResidueTypeValue,
        equilibration_mode=equilibrationModeValue,
        spawning=spawningValue,
        log_file=logFileValue,
        rescoring=rescoringValue,
        epsilon=epsilonValue,
        ligand_equilibration_cst=ligandEquilibrationCstValue,

        
        # Implement all the variables...
    )

    from utils import launchCalculationAction

    launchCalculationAction(
        block,
        jobs,
        program="pele",
        uploadFolders=[
            "pele",
            poses_folder,
        ],
    )


def peleFinalAction(block: SlurmBlock):
    print("Pele finished")

    from utils import downloadResultsAction

    downloadResultsAction(block)

    peleFolderName = block.variables.get("pele_folder_name", "pele")

    block.setOutput("pele_output_folder", peleFolderName)


from utils import BSC_JOB_VARIABLES

blockVariables = BSC_JOB_VARIABLES + [
    peleFolderNameVariable,
    selectionsListVariable,
    boxCentersVariable,
    constraintsVariable,
    ligandIndexVariable,
    peleStepsVariable,
    peleDebugVariable,
    peleIterationsVariable,
    equilibrationStepsVariable,
    ligandEnergyGroupsVariable,
    peleSeparatorVariable,
    usePeleffyVariable,
    useSrunVariable,
    energyByResidueVariable,
    ebrNewFlagVariable,
    ninetyDegreesVersionVariable,
    analysisVariable,
    energyByResidueTypeVariable,
    peptideVariable,
    equilibrationModeVariable,
    spawningVariable,
    continuationVariable,
    equilibrationVariable,
    skipModelsVariable,
    skipLigandsVariable,
    extendIterationsVariable,
    onlyModelsVariable,
    onlyLigandsVariable,
    onlyCombinationsVariable,
    ligandTemplateVariable,
    seedVariable,
    logFileVariable,
    rescoringVariable,
    epsilonVariable
]

peleBlock = SlurmBlock(
    name="PELE",
    description="Run PELE",
    initialAction=peleAction,
    finalAction=peleFinalAction,
    inputGroups=[folderInputGroup, glideOutputGroup],
    variables=blockVariables,
    outputs=[peleOutputFolderOutput],
)
