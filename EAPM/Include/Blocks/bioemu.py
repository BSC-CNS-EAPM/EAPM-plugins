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
the BioEmu conda environment and exports COLABFOLD_DIR, so the program name is
passed through the shared EAPM launcher.

The GPU request is exposed as the "GPUs" Slurm variable (``--gres gpu:N``, only
written on the ``acc_*`` partitions) and the CPUs that go with them are set with
the shared "CPUs per task" variable, e.g. 1 GPU + 20 CPUs per task.
"""

from utils import BSC_JOB_VARIABLES, gpusVariable

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
clusterEnvVariable = PluginVariable(
    name="Cluster BioEmu environment",
    id="cluster_env",
    description="Conda environment activated by the job on the cluster, as a path or "
    "a name. Leave it empty to use the one bsc_calculations hardcodes for BioEmu "
    "(/gpfs/projects/bsc72/conda_envs/bioemu).",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/bsc72/conda_envs/bioemu",
)
colabfoldBinVariable = PluginVariable(
    name="ColabFold bin folder",
    id="colabfold_bin",
    description="Folder holding 'colabfold_search', prepended to the job's PATH. "
    "Adapt it to the ColabFold installation that goes with your environment.",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/bsc72/conda_envs/bioemu/colabfold/localcolabfold/colabfold-conda/bin",
)
colabfoldDirVariable = PluginVariable(
    name="COLABFOLD_DIR",
    id="colabfold_dir",
    description="Exported as COLABFOLD_DIR when set. Leave it empty if your ColabFold "
    "installation does not need it.",
    type=VariableTypes.STRING,
    defaultValue="",
)
hfCacheVariable = PluginVariable(
    name="HuggingFace cache",
    id="hf_cache",
    description="Pre-populated HuggingFace cache holding the BioEmu checkpoints. "
    "Exported as HF_HOME (with HF_HUB_CACHE=<cache>/hub) together with the offline "
    "flags, so the compute node never reaches out for the model weights. Empty "
    "disables the exports, which needs the node to have internet access.",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/bsc72/bioemu_hf_cache",
)
clusterModulesVariable = PluginVariable(
    name="Cluster modules",
    id="cluster_modules",
    description="Comma-separated modules loaded by the job before activating the "
    "environment.",
    type=VariableTypes.STRING,
    defaultValue="bsc/1.0, anaconda, intel/2025.1",
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
    description="Build the MSA on the cluster with colabfold_search against a local "
    "database. Leave it off only if the compute nodes have internet access: BioEmu "
    "otherwise queries the ColabFold web API to align the raw sequence.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)
mmseqsPathVariable = PluginVariable(
    name="MMseqs2 binary",
    id="mmseqs_path",
    description="Path to the MMseqs2 executable used by colabfold_search "
    "(only used when 'MSA calculation' is enabled).",
    type=VariableTypes.STRING,
    defaultValue="/apps/ACC/MMSEQS2/17-b804f/GCC/OPENMPI/bin/mmseqs",
)
msaDatabaseVariable = PluginVariable(
    name="ColabFold database",
    id="msa_database",
    description="Path to the local ColabFold/MMseqs2 database searched for the MSA "
    "(only used when 'MSA calculation' is enabled).",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/apps/MN5/ACC/COLABFOLD/SRC/database/FULL",
)
clusterMsaFolderVariable = PluginVariable(
    name="Persistent MSA folder",
    id="cluster_msa_folder",
    description="Absolute folder on the cluster where the alignments are kept, one "
    "subfolder per model. Living outside the folder uploaded on each launch, an "
    "alignment found there is reused and colabfold_search is skipped. Empty keeps "
    "them inside the run folder, so every launch searches again.",
    type=VariableTypes.STRING,
    defaultValue="",
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
        raise Exception("The JSON sequences file must map sequence names to sequences.")

    fasta_path = os.path.join(os.getcwd(), "bioemu_input.fasta")
    with open(fasta_path, "w") as of:
        for name, sequence in data.items():
            of.write(f">{name}\n{sequence}\n")
    print(f"Converted the JSON sequences into {fasta_path}")
    return fasta_path


def _use_local_msa(jobs, folder_name, mmseqs_path, database_path, msa_root=None):
    """
    Prepend a local ``colabfold_search`` to each job and feed it its MSA.

    Without this BioEmu receives the raw sequence and aligns it through the
    ColabFold web API, which the compute nodes cannot reach. ``setUpBioEmu`` has
    an ``msa_calculation`` flag meant to do this, but it writes the search and
    the ``RUN_SAMPLES`` assignment on the same line and then reads the alignment
    from a path the search never produces, so the step is built here instead.

    ``msa_root`` keeps the alignments outside the folder uploaded on each launch,
    so a rerun reuses them instead of searching again.
    """
    import os
    import re

    msa_jobs = []
    for job in jobs:
        model_folder = re.search(r"--output_dir (\S+)", job)
        if model_folder is None:
            raise Exception("Could not tell which model a BioEmu job belongs to.")

        model = os.path.basename(model_folder.group(1).rstrip("/"))
        fasta_file = os.path.join(folder_name, "fastas", model + ".fasta")
        msa_folder = os.path.join(msa_root or os.path.join(folder_name, "msas"), model)

        # colabfold_search names the alignment after the FASTA header, so pick up
        # whatever .a3m it wrote rather than guessing the name.
        search = (
            f"mkdir -p {msa_folder}\n"
            f"if ! ls {msa_folder}/*.a3m > /dev/null 2>&1; then\n"
            f"colabfold_search --mmseqs {mmseqs_path} "
            f"{fasta_file} {database_path} {msa_folder}\n"
            f"fi\n"
            f"MSA=$(ls {msa_folder}/*.a3m | head -n 1)\n"
        )

        msa_jobs.append(search + re.sub(r"--sequence \S+ ", '--sequence "$MSA" ', job, count=1))

    return msa_jobs


def _use_env_python(jobs, cluster_env):
    """
    Call the environment's interpreter by path instead of relying on ``python``.

    The ColabFold installation ships its own interpreter and has to sit at the
    front of PATH for ``colabfold_search``, which leaves a bare ``python``
    pointing at ColabFold rather than at the activated BioEmu environment.
    """
    import os
    import re

    env_python = os.path.join(cluster_env, "bin", "python")

    return [re.sub(r"(?<![\w./-])python(?= )", env_python, job) for job in jobs]


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
        # The MSA step is added below, see _use_local_msa.
        msa_calculation=False,
        bioemu_env=block.variables.get(bioemuEnvVariable.id) or None,
        conda_sh=block.variables.get(condaShVariable.id) or "~/miniconda3/etc/profile.d/conda.sh",
    )

    if not jobs:
        raise Exception(
            "No BioEmu jobs were generated. If 'Skip finished' is enabled the "
            "sampling may already be complete."
        )

    if block.variables.get(msaCalculationVariable.id, False):
        print("Adding the local MSA calculation to the jobs...")
        msa_root = (block.variables.get(clusterMsaFolderVariable.id) or "").strip()
        if msa_root:
            print(f"Reusing any alignment already present in '{msa_root}'")
        jobs = _use_local_msa(
            jobs,
            folder_name,
            block.variables.get(mmseqsPathVariable.id)
            or "/apps/ACC/MMSEQS2/17-b804f/GCC/OPENMPI/bin/mmseqs",
            block.variables.get(msaDatabaseVariable.id)
            or "/gpfs/apps/MN5/ACC/COLABFOLD/SRC/database/FULL",
            msa_root=msa_root or None,
        )

    print(f"Generated {len(jobs)} BioEmu job(s).")

    # bsc_calculations pins the BioEmu environment whenever it is told the program
    # is 'bioemu', so a custom one means setting up the environment here instead.
    cluster_env = (block.variables.get(clusterEnvVariable.id) or "").strip()

    if not cluster_env:
        launchCalculationAction(block, jobs, "bioemu", [folder_name])
        return

    jobs = _use_env_python(jobs, cluster_env)

    exports = []

    colabfold_bin = (block.variables.get(colabfoldBinVariable.id) or "").strip()
    if colabfold_bin:
        exports.append(f"PATH={colabfold_bin}:$PATH")

    colabfold_dir = (block.variables.get(colabfoldDirVariable.id) or "").strip()
    if colabfold_dir:
        exports.append(f"COLABFOLD_DIR={colabfold_dir}")

    # BioEmu pulls its checkpoints from the HuggingFace hub on first use. Pointing it
    # at a populated cache and forcing offline mode keeps that off the network, which
    # the compute nodes cannot use anyway.
    hf_cache = (block.variables.get(hfCacheVariable.id) or "").strip()
    if hf_cache:
        exports += [
            f"HF_HOME={hf_cache}",
            f"HF_HUB_CACHE={hf_cache}/hub",
            "HF_HUB_OFFLINE=1",
            "TRANSFORMERS_OFFLINE=1",
        ]

    modules = [
        module.strip()
        for module in (block.variables.get(clusterModulesVariable.id) or "").split(",")
        if module.strip()
    ]

    print(f"Using the cluster environment '{cluster_env}'")
    if hf_cache:
        print(f"Reading the BioEmu checkpoints from the cache at '{hf_cache}'")

    launchCalculationAction(
        block,
        jobs,
        None,
        [folder_name],
        condaEnv=cluster_env,
        modules=modules or None,
        exports=exports,
    )


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
        gpusVariable,
        output,
        numSamplesVariable,
        batchSizeVariable,
        bioemuEnvVariable,
        clusterEnvVariable,
        clusterModulesVariable,
        colabfoldBinVariable,
        colabfoldDirVariable,
        hfCacheVariable,
        condaShVariable,
        gpuLocalVariable,
        filterSamplesVariable,
        msaCalculationVariable,
        mmseqsPathVariable,
        msaDatabaseVariable,
        clusterMsaFolderVariable,
        skipFinishedVariable,
        removeExistingResults,
    ],
    inputs=[sequencesFile],
    outputs=[bioemuFolderVariable],
)
