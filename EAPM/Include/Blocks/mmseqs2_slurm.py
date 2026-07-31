"""
Run MMseqs2 sequence clustering on an HPC cluster via SLURM.

A runner script (SBATCH headers + the ``mmseqs easy-cluster`` command) is
generated in the flow working directory. For a remote, the input FASTA and
runner are transferred to a per-run sandbox on the cluster and the job is
submitted there; the final action downloads the sandbox back and builds the
``{representative: [members]}`` clustering map. For the local remote the runner
is submitted directly with no transfer.

This is the SLURM counterpart of the local ``MMseqs2 Clustering`` block; it runs
the exact same ``mmseqs easy-cluster`` call.
"""

import os
import posixpath
import shlex

from HorusAPI import PluginVariable, SlurmBlock, VariableTypes

# ==========================#
# Input
# ==========================#
sequencesFile = PluginVariable(
    id="sequences_file",
    name="Sequences file",
    description="Sequences to cluster. Either a FASTA file or a JSON file mapping "
    "sequence name to sequence (e.g. upo_dataset_sequences.json).",
    type=VariableTypes.FILE,
    allowedValues=["fasta", "fa", "faa", "json"],
)

# ==========================#
# MMseqs2 parameters
# ==========================#
minSeqIdVar = PluginVariable(
    id="min_seq_id",
    name="Minimum sequence identity",
    description="List matches above this sequence identity (--min-seq-id, 0.0-1.0).",
    type=VariableTypes.FLOAT,
    defaultValue=0.5,
)
coverageVar = PluginVariable(
    id="coverage",
    name="Coverage",
    description="Minimum alignment coverage (-c, 0.0-1.0).",
    type=VariableTypes.FLOAT,
    defaultValue=0.8,
)
covModeVar = PluginVariable(
    id="cov_mode",
    name="Coverage mode",
    description="MMseqs2 coverage mode (--cov-mode). 0: bidirectional, "
    "1: target coverage, 2: query coverage.",
    type=VariableTypes.INTEGER,
    defaultValue=1,
)

# ==========================#
# SLURM parameters
# ==========================#
preambleTextVar = PluginVariable(
    id="preamble",
    name="Preamble",
    description=(
        "Shell lines prepended to the runner before MMseqs2 runs, to put 'mmseqs' "
        "on the PATH (e.g. 'module load mmseqs2'). One command per line. "
        "Combined with the Preamble file below if both are set."
    ),
    type=VariableTypes.TEXT_AREA,
    defaultValue=None,
)
preambleFileVar = PluginVariable(
    id="preamble_file",
    name="Preamble file",
    description=(
        "Shell snippet prepended to the runner (module loads, conda activate, ...). "
        "Typically puts 'mmseqs' on the PATH (e.g. 'module load mmseqs2'). "
        "Leave empty for none."
    ),
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["sh", "txt"],
)
ntasksVar = PluginVariable(
    id="ntasks",
    name="ntasks",
    description="#SBATCH --ntasks value.",
    type=VariableTypes.INTEGER,
    defaultValue=1,
)
cpusPerTaskVar = PluginVariable(
    id="cpus_per_task",
    name="cpus-per-task",
    description="#SBATCH --cpus-per-task value. Also passed to MMseqs2 as --threads.",
    type=VariableTypes.INTEGER,
    defaultValue=4,
)
timeVar = PluginVariable(
    id="time",
    name="Time limit",
    description="#SBATCH --time value (e.g. 00-01:00:00).",
    type=VariableTypes.STRING,
    defaultValue="00-01:00:00",
)
qosVar = PluginVariable(
    id="qos",
    name="QOS",
    description="#SBATCH --qos value. Leave empty for none.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
partitionVar = PluginVariable(
    id="partition",
    name="Partition",
    description="#SBATCH --partition value. Leave empty for none.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
accountVar = PluginVariable(
    id="account",
    name="Account",
    description="#SBATCH --account value. Leave empty for none.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
jobNameVar = PluginVariable(
    id="job_name",
    name="Job name",
    description="SLURM job name. Leave empty to derive from the input file name.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
mmseqsCommandVar = PluginVariable(
    id="mmseqs_command",
    name="MMseqs2 command",
    description="Command used to invoke MMseqs2 on the cluster (after the preamble "
    "has loaded it). Usually just 'mmseqs'.",
    type=VariableTypes.STRING,
    defaultValue="mmseqs",
)

# ==========================#
# Output
# ==========================#
clustersFile = PluginVariable(
    id="clusters_file",
    name="Clusters file",
    description="JSON file mapping each cluster representative to the list of its members.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
representativesFile = PluginVariable(
    id="representatives_file",
    name="Representatives FASTA",
    description="FASTA file with one representative sequence per cluster.",
    type=VariableTypes.FILE,
    allowedValues=["fasta"],
)
resultsFolder = PluginVariable(
    id="results_folder",
    name="MMseqs2 results folder",
    description="Folder holding the raw MMseqs2 output (cluster TSV, representative "
    "and all-sequences FASTA).",
    type=VariableTypes.FOLDER,
)

# MMseqs2 easy-cluster output prefix used inside the sandbox
CLUSTER_PREFIX = "clusterRes"


def _ensure_local_fasta(sequences_path: str, base: str) -> str:
    """Return a path to a FASTA file, converting a JSON ``{name: seq}`` file if needed."""
    import json

    if not sequences_path.lower().endswith(".json"):
        return os.path.abspath(sequences_path)

    with open(sequences_path) as jf:
        data = json.load(jf)
    if not isinstance(data, dict):
        raise ValueError("The JSON sequences file must map sequence names to sequences.")

    fasta_path = os.path.join(base, "mmseqs_input.fasta")
    with open(fasta_path, "w") as of:
        for name, seq in data.items():
            of.write(f">{name}\n{seq}\n")
    return fasta_path


def _sbatch_header(block, job_name: str) -> str:
    lines = ["#!/bin/bash", f"#SBATCH --job-name={job_name}"]

    ntasks = block.variables.get("ntasks")
    if ntasks not in (None, ""):
        lines.append(f"#SBATCH --ntasks={int(ntasks)}")

    cpus = block.variables.get("cpus_per_task")
    if cpus not in (None, ""):
        lines.append(f"#SBATCH --cpus-per-task={int(cpus)}")

    time_limit = block.variables.get("time")
    if time_limit:
        lines.append(f"#SBATCH --time={time_limit}")

    qos = block.variables.get("qos")
    if qos:
        lines.append(f"#SBATCH --qos={qos}")

    partition = block.variables.get("partition")
    if partition:
        lines.append(f"#SBATCH --partition={partition}")

    account = block.variables.get("account")
    if account:
        lines.append(f"#SBATCH --account={account}")

    lines.append("#SBATCH --output=mmseqs_%j.out")
    lines.append("#SBATCH --error=mmseqs_%j.err")
    return "\n".join(lines)


def _preamble_text(block) -> str:
    """Combine the inline preamble and the preamble file (in that order)."""
    parts = []

    inline = block.variables.get("preamble")
    if inline and str(inline).strip():
        parts.append(str(inline).rstrip())

    preamble_file = block.variables.get("preamble_file")
    if preamble_file:
        with open(os.path.expanduser(str(preamble_file)), "r", encoding="utf-8") as fh:
            parts.append(fh.read().rstrip())

    if not parts:
        return ""
    return "\n".join(parts) + "\n"


def _build_command(block, input_basename: str) -> str:
    """Build the ``mmseqs easy-cluster`` command run inside the sandbox."""
    mmseqs_cmd = block.variables.get("mmseqs_command") or "mmseqs"

    parts = [
        shlex.quote(mmseqs_cmd),
        "easy-cluster",
        shlex.quote(input_basename),
        CLUSTER_PREFIX,
        "tmp",
        "--min-seq-id",
        str(block.variables.get("min_seq_id", 0.5)),
        "-c",
        str(block.variables.get("coverage", 0.8)),
        "--cov-mode",
        str(block.variables.get("cov_mode", 1)),
    ]

    cpus = block.variables.get("cpus_per_task")
    if cpus not in (None, ""):
        parts += ["--threads", str(int(cpus))]

    return " ".join(parts)


def submitMMseqsSlurm(block: SlurmBlock):
    sequences_path = block.inputs.get("sequences_file")
    if not sequences_path or not os.path.isfile(sequences_path):
        raise Exception("A valid sequences file (FASTA or JSON) must be provided.")

    base = os.getcwd()
    local_fasta = _ensure_local_fasta(sequences_path, base)
    input_basename = os.path.basename(local_fasta)

    stem = os.path.splitext(input_basename)[0]
    job_name = block.variables.get("job_name") or f"mmseqs_{stem}"

    header = _sbatch_header(block, job_name)
    preamble = _preamble_text(block)
    command = _build_command(block, input_basename)

    # ---------------------------------------------------------------- #
    # Local remote: run in place, no transfer.
    # ---------------------------------------------------------------- #
    if block.remote.isLocal:
        input_in_base = os.path.join(base, input_basename)
        if os.path.abspath(local_fasta) != input_in_base:
            block.remote.command(
                f"cp {shlex.quote(os.path.abspath(local_fasta))} {shlex.quote(input_in_base)}"
            )

        runner_text = f"{header}\n\n{preamble}cd {shlex.quote(base)}\n{command}\n"
        runner = os.path.join(base, f"{job_name}_runner.sh")
        with open(runner, "w", encoding="utf-8") as fh:
            fh.write(runner_text)

        block.extraData["results_folder"] = base
        block.extraData["local_results_folder"] = base

        print(f"Submitting MMseqs2 SLURM job: {runner}")
        job_id = block.remote.submitJob(runner)
        print(f"Submitted MMseqs2 SLURM job: {job_id}")
        return

    # ---------------------------------------------------------------- #
    # Remote: transfer the FASTA + runner to a sandbox and submit.
    # ---------------------------------------------------------------- #
    remote_sandbox = posixpath.join(block.remote.workDir, job_name)
    block.remote.command(f"mkdir -p {shlex.quote(remote_sandbox)}")

    print(f"Transferring sequences to '{remote_sandbox}'...")
    block.remote.sendData(local_fasta, remote_sandbox)

    runner_text = f"{header}\n\n{preamble}cd {shlex.quote(remote_sandbox)}\n{command}\n"
    runner = os.path.join(base, f"{job_name}_runner.sh")
    with open(runner, "w", encoding="utf-8") as fh:
        fh.write(runner_text)

    remote_runner = block.remote.sendData(runner, remote_sandbox)
    print(f"Transferred runner to: {remote_runner}")

    local_sandbox = os.path.join(base, job_name)
    block.extraData["results_folder"] = remote_sandbox
    block.extraData["local_results_folder"] = local_sandbox

    print(f"Submitting MMseqs2 SLURM job: {remote_runner}")
    job_id = block.remote.submitJob(remote_runner)
    print(f"Submitted MMseqs2 SLURM job: {job_id}")


def collectMMseqsSlurm(block: SlurmBlock):
    import json
    import shutil

    remote_results = block.extraData.get("results_folder")
    local_results = block.extraData.get("local_results_folder", remote_results)

    if not block.remote.isLocal:
        print(f"Downloading results from '{remote_results}'...")
        local_results = block.remote.getData(remote_results, local_results)

    base = os.getcwd()
    cluster_tsv = os.path.join(local_results, CLUSTER_PREFIX + "_cluster.tsv")
    if not os.path.isfile(cluster_tsv):
        raise Exception(
            f"MMseqs2 clustering result not found ({cluster_tsv}). "
            "Check the mmseqs_*.err file in the results folder."
        )

    clusters: dict = {}
    with open(cluster_tsv) as cr:
        for line in cr:
            if not line.strip():
                continue
            representative, member = line.split()
            clusters.setdefault(representative, []).append(member)

    clusters_output = os.path.join(base, "mmseqs_clusters.json")
    with open(clusters_output, "w") as jf:
        json.dump(clusters, jf, indent=2)

    representatives_output = os.path.join(base, "mmseqs_representatives.fasta")
    rep_seq_fasta = os.path.join(local_results, CLUSTER_PREFIX + "_rep_seq.fasta")
    if os.path.isfile(rep_seq_fasta):
        shutil.copyfile(rep_seq_fasta, representatives_output)

    print(f"MMseqs2 SLURM run finished. Found {len(clusters)} clusters. "
          f"Results in: {local_results}")

    block.setOutput(clustersFile.id, clusters_output)
    if os.path.isfile(representatives_output):
        block.setOutput(representativesFile.id, representatives_output)
    block.setOutput(resultsFolder.id, local_results)


mmseqsClusterSlurmBlock = SlurmBlock(
    name="MMseqs2 Clustering (SLURM)",
    id="mmseqs2_cluster_slurm",
    description="Cluster protein sequences by identity with MMseqs2 (easy-cluster), "
    "submitted as a SLURM job.",
    initialAction=submitMMseqsSlurm,
    finalAction=collectMMseqsSlurm,
    inputs=[sequencesFile],
    variables=[
        minSeqIdVar,
        coverageVar,
        covModeVar,
        preambleTextVar,
        preambleFileVar,
        ntasksVar,
        cpusPerTaskVar,
        timeVar,
        qosVar,
        partitionVar,
        accountVar,
        jobNameVar,
        mmseqsCommandVar,
    ],
    outputs=[clustersFile, representativesFile, resultsFolder],
    category="Clustering & Selection",
)
