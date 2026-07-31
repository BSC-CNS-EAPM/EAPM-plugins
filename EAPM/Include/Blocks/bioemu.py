"""
Module containing the BioEmu sampling block for the EAPM plugin.

Reproduces the UPO notebook's BioEmu setup and submission (cells 38-40)::

    jobs = selected_sequences.setUpBioEmu('bioemu_sampling', num_samples=10000,
                                          batch_size_100=200, bioemu_env='bioemu',
                                          gpu_local=False, skip_finished=True)
    bsc_calculations.mn5.jobArrays(jobs, program='bioemu', gpus=1, ...)

BioEmu runs off ``prepare_proteins.sequenceModels`` (sequences, not structures),
so the block takes a FASTA of the selected sequences.

``bsc_calculations`` special-cases ``program='bioemu'`` for MareNostrum: it loads
the BioEmu conda environment, exports COLABFOLD_DIR and requests a GPU
(``--gres gpu:1`` by default), so nothing extra is needed here beyond passing the
program name through the shared EAPM launcher.
"""

from utils import BSC_JOB_VARIABLES

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
sequencesFile = PluginVariable(
    name="Sequences file",
    id="sequences_file",
    description="Sequences to sample. A FASTA file (e.g. the combined selection) or "
    "a JSON {name: sequence} mapping.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta", "fa", "faa", "json"],
)

# ==========================#
# Variables
# ==========================#
output = PluginVariable(
    name="BioEmu folder",
    id="folder_name",
    description="Name of the folder where the BioEmu sampling will be stored.",
    type=VariableTypes.STRING,
    defaultValue="bioemu_sampling",
)
numSamplesVariable = PluginVariable(
    name="Number of samples",
    id="num_samples",
    description="Number of conformations sampled per model.",
    type=VariableTypes.INTEGER,
    defaultValue=10000,
)
batchSizeVariable = PluginVariable(
    name="Batch size (per 100 residues)",
    id="batch_size_100",
    description="BioEmu batch size normalised to 100 residues. The notebook uses 200; "
    "the prepare_proteins default is 20. Lower it if you hit GPU memory limits.",
    type=VariableTypes.INTEGER,
    defaultValue=200,
)
bioemuEnvVariable = PluginVariable(
    name="BioEmu environment",
    id="bioemu_env",
    description="Name of the conda environment holding BioEmu.",
    type=VariableTypes.STRING,
    defaultValue="bioemu",
)
condaShVariable = PluginVariable(
    name="conda.sh path",
    id="conda_sh",
    description="Path to the conda profile script used to activate the environment.",
    type=VariableTypes.STRING,
    defaultValue="~/miniconda3/etc/profile.d/conda.sh",
)
gpuLocalVariable = PluginVariable(
    name="Local GPU",
    id="gpu_local",
    description="Build the commands for a local GPU run instead of a cluster job.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
filterSamplesVariable = PluginVariable(
    name="Filter samples",
    id="filter_samples",
    description="Apply BioEmu's sample filtering.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)
msaCalculationVariable = PluginVariable(
    name="MSA calculation",
    id="msa_calculation",
    description="Compute MSAs for the sequences (requires the MMseqs2/ColabFold "
    "databases configured on the cluster).",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
skipFinishedVariable = PluginVariable(
    name="Skip finished",
    id="skip_finished",
    description="Do not regenerate jobs whose sampling already completed.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)
removeExistingResults = PluginVariable(
    name="Remove existing results",
    id="remove_existing_results",
    description="Remove the BioEmu folder if it already exists.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)

# ==========================#
# Variable outputs
# ==========================#
bioemuFolderVariable = PluginVariable(
    id="bioemu_folder",
    name="BioEmu folder",
    description="Folder holding the BioEmu sampling output (feed this into the "
    "BioEmu analysis blocks).",
    type=VariableTypes.FOLDER,
)


def _ensure_fasta(sequences_path: str) -> str:
    """Return a FASTA path, converting a JSON ``{name: sequence}`` file if needed."""
    import json
    import os

    if not sequences_path.lower().endswith(".json"):
        return sequences_path

    with open(sequences_path) as jf:
        data = json.load(jf)
    if not isinstance(data, dict):
        raise Exception(
            "The JSON sequences file must map sequence names to sequences."
        )

    fasta_path = os.path.join(os.getcwd(), "bioemu_input.fasta")
    with open(fasta_path, "w") as of:
        for name, sequence in data.items():
            of.write(f">{name}\n{sequence}\n")
    print(f"Converted the JSON sequences into {fasta_path}")
    return fasta_path


def initial_bioemu(block: SlurmBlock):
    """
    Set up the BioEmu sampling jobs and send them to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """
    # pylint: disable=import-outside-toplevel
    import os

    import prepare_proteins
    from utils import launchCalculationAction

    # pylint: enable=import-outside-toplevel

    sequences_file = block.inputs.get(sequencesFile.id, None)
    if not sequences_file or sequences_file == "None":
        raise Exception("No sequences file provided.")
    if not os.path.isfile(sequences_file):
        raise Exception(f"The sequences file '{sequences_file}' does not exist.")

    sequences_file = _ensure_fasta(sequences_file)

    folder_name = block.variables.get(output.id, "bioemu_sampling")
    remove_existing = block.variables.get(removeExistingResults.id, False)

    if remove_existing and os.path.exists(folder_name):
        os.system("rm -rf " + folder_name)

    if not remove_existing and os.path.exists(folder_name):
        raise Exception(
            f"The folder {folder_name} already exists. "
            "Please, choose another name or remove it with the remove existing folder option."
        )

    block.extraData["folder_name"] = folder_name

    print("Loading sequences...")
    sequences = prepare_proteins.sequenceModels(sequences_file)

    print("Setting up BioEmu...")
    jobs = sequences.setUpBioEmu(
        folder_name,
        num_samples=block.variables.get(numSamplesVariable.id, 10000),
        batch_size_100=block.variables.get(batchSizeVariable.id, 200),
        gpu_local=block.variables.get(gpuLocalVariable.id, False),
        skip_finished=block.variables.get(skipFinishedVariable.id, True),
        filter_samples=block.variables.get(filterSamplesVariable.id, True),
        msa_calculation=block.variables.get(msaCalculationVariable.id, False),
        bioemu_env=block.variables.get(bioemuEnvVariable.id) or None,
        conda_sh=block.variables.get(condaShVariable.id)
        or "~/miniconda3/etc/profile.d/conda.sh",
    )

    if not jobs:
        raise Exception(
            "No BioEmu jobs were generated. If 'Skip finished' is enabled the "
            "sampling may already be complete."
        )

    print(f"Generated {len(jobs)} BioEmu job(s).")

    launchCalculationAction(block, jobs, "bioemu", [folder_name])


def final_bioemu(block: SlurmBlock):
    """
    Download the BioEmu sampling results.

    Args:
        block (SlurmBlock): The SlurmBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import os

    from utils import downloadResultsAction

    # pylint: enable=import-outside-toplevel

    downloaded_path = downloadResultsAction(block)

    results_folder = block.extraData["folder_name"]
    bioemu_folder = os.path.join(downloaded_path, results_folder)

    print(f"BioEmu sampling finished. Results in: {bioemu_folder}")

    block.setOutput(bioemuFolderVariable.id, bioemu_folder)


bioEmuBlock = SlurmBlock(
    category="BioEmu",
    name="BioEmu Sampling",
    id="bioemu_sampling",
    description="Sample conformational ensembles with BioEmu "
    "(for MareNostrum, Nord3 clusters or local).",
    initialAction=initial_bioemu,
    finalAction=final_bioemu,
    variables=BSC_JOB_VARIABLES
    + [
        output,
        numSamplesVariable,
        batchSizeVariable,
        bioemuEnvVariable,
        condaShVariable,
        gpuLocalVariable,
        filterSamplesVariable,
        msaCalculationVariable,
        skipFinishedVariable,
        removeExistingResults,
    ],
    inputs=[sequencesFile],
    outputs=[bioemuFolderVariable],
)
