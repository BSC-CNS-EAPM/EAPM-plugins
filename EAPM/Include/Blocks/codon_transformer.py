"""
Module containing the CodonTransformer block for the EAPM plugin.

Reproduces the UPO notebook's codon optimisation (cells 71-74)::

    codon_sequences = bioprospecting.expression.codonTransformer.computeCodonTransformerSequences(
        selected_df.to_dict()['Sequence'],
        organism="Escherichia coli general",
        device="cuda",
        codon_transformer_conda_env="codontransformer")

Given the protein sequences of the selected models, CodonTransformer returns a
DNA sequence optimised for expression in the requested host.

CodonTransformer runs in its own conda environment (it pulls in torch and the
transformer weights), which is why the environment name is a parameter rather
than a plugin dependency.

The environment is resolved to its ``bin`` directory and prepended to PATH so the
subprocess picks up the right interpreter. This avoids ``conda activate``, which is
unavailable on hosts that use micromamba.
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
sequencesFile = PluginVariable(
    id="sequences_file",
    name="Protein sequences",
    description="Protein sequences to codon-optimise. A FASTA file or a JSON "
    "{name: sequence} mapping.",
    type=VariableTypes.FILE,
    allowedValues=["fasta", "fa", "faa", "json"],
)
modelsFile = PluginVariable(
    id="models_file",
    name="Models subset (optional)",
    description="Optional JSON list of model names (e.g. the Pareto selection) used "
    "to restrict the optimisation to those sequences.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
organismVariable = PluginVariable(
    id="organism",
    name="Host organism",
    description="Expression host whose codon usage is targeted.",
    type=VariableTypes.STRING,
    defaultValue="Escherichia coli general",
)
deviceVariable = PluginVariable(
    id="device",
    name="Device",
    description="Torch device used to run the model. Use 'cuda' when a GPU is "
    "available; 'cpu' works but is slower.",
    type=VariableTypes.STRING,
    defaultValue="cpu",
    allowedValues=["cpu", "cuda"],
)
condaEnvVariable = PluginVariable(
    id="codon_transformer_conda_env",
    name="CodonTransformer environment",
    description="Environment holding CodonTransformer: an env name (looked up in the "
    "conda/micromamba roots), a path to the environment, or a path to its python. "
    "Leave empty to use the Horus backend's own interpreter.",
    type=VariableTypes.STRING,
    defaultValue="codontransformer",
)
startMethionineVariable = PluginVariable(
    id="add_start_methionine",
    name="Add start methionine",
    description="Prepend a start methionine to sequences that lack one.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)

# ==========================#
# Variable outputs
# ==========================#
dnaFastaFile = PluginVariable(
    id="dna_fasta_file",
    name="Optimised DNA (FASTA)",
    description="FASTA with the codon-optimised DNA sequences.",
    type=VariableTypes.FILE,
    allowedValues=["fasta"],
)
dnaJsonFile = PluginVariable(
    id="dna_json_file",
    name="Optimised DNA (JSON)",
    description="JSON mapping each model to its optimised DNA sequence.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
tableFile = PluginVariable(
    id="table_file",
    name="Expression table",
    description="CSV with the protein sequence, the optimised DNA and its length "
    "for each model.",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)


def _read_sequences(sequences_path: str) -> dict:
    """Read sequences from a FASTA or JSON file into a ``{name: sequence}`` dict."""
    import json

    if sequences_path.lower().endswith(".json"):
        with open(sequences_path) as jf:
            data = json.load(jf)
        if not isinstance(data, dict):
            raise ValueError(
                "The JSON sequences file must map sequence names to sequences."
            )
        return {str(name): str(seq) for name, seq in data.items()}

    sequences: dict = {}
    name = None
    chunks: list = []
    with open(sequences_path) as ff:
        for line in ff:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    if not sequences:
        raise ValueError("No sequences found in the sequences file.")
    return sequences


def _resolve_environment_bin(env_value: str):
    """
    Resolve the CodonTransformer environment to the ``bin`` directory holding its python.

    Accepts an environment name (searched for in the usual micromamba/conda roots), a
    path to the environment itself, or a path to its python interpreter. Returns None
    when no environment is given, meaning the backend's own interpreter is used.
    """
    # pylint: disable=import-outside-toplevel
    import os

    # pylint: enable=import-outside-toplevel

    if not env_value or not env_value.strip():
        return None
    env_value = env_value.strip()

    if os.sep in env_value:
        candidates = [os.path.abspath(os.path.expanduser(env_value))]
    else:
        roots = [
            os.environ.get("MAMBA_ROOT_PREFIX"),
            os.path.expanduser("~/micromamba"),
            os.path.expanduser("~/miniforge3"),
            os.path.expanduser("~/miniconda3"),
            os.path.expanduser("~/anaconda3"),
        ]
        conda_exe = os.environ.get("CONDA_EXE")
        if conda_exe:
            roots.append(os.path.dirname(os.path.dirname(conda_exe)))
        candidates = list(
            dict.fromkeys(
                os.path.join(root, "envs", env_value) for root in roots if root
            )
        )

    for candidate in candidates:
        # A path to the interpreter itself.
        if os.path.isfile(candidate):
            return os.path.dirname(candidate)
        # An environment prefix.
        if os.path.isfile(os.path.join(candidate, "bin", "python")):
            return os.path.join(candidate, "bin")

    raise ValueError(
        f"Could not find the CodonTransformer environment '{env_value}'. Give an "
        "existing environment name, the path to the environment, or the path to its "
        f"python. Looked in: {', '.join(candidates)}"
    )


def _assert_codon_transformer(env_bin, env_value: str):
    """Fail early, with a readable message, if CodonTransformer is not importable."""
    # pylint: disable=import-outside-toplevel
    import os
    import subprocess
    import sys

    # pylint: enable=import-outside-toplevel

    python = os.path.join(env_bin, "python") if env_bin else sys.executable
    probe = subprocess.run(
        [python, "-c", "import CodonTransformer"],
        capture_output=True,
        text=True,
        check=False,
    )
    if probe.returncode != 0:
        where = f"environment '{env_value}'" if env_value else "the Horus backend"
        raise ValueError(
            f"CodonTransformer is not installed in {where} ({python}). Install it there "
            f"with 'pip install CodonTransformer'. Error: {probe.stderr.strip()}"
        )
    return python


def compute_codon_sequences(block: PluginBlock):
    """
    Codon-optimise the input protein sequences for the requested host.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import csv
    import json
    import os

    # pylint: enable=import-outside-toplevel

    if block.remote.name != "Local":
        raise ValueError("This block only works on the local machine.")

    sequences_path = block.inputs.get(sequencesFile.id, None)
    if not sequences_path or not os.path.isfile(sequences_path):
        raise ValueError("A valid protein sequences file (FASTA or JSON) must be provided.")

    sequences = _read_sequences(sequences_path)

    models_path = block.inputs.get(modelsFile.id, None)
    if models_path and os.path.isfile(models_path):
        with open(models_path) as jf:
            wanted = json.load(jf)
        if isinstance(wanted, dict):
            wanted = list(wanted)
        wanted = [str(m) for m in wanted]
        missing = [m for m in wanted if m not in sequences]
        if missing:
            print(f"Warning: {len(missing)} requested model(s) are not in the "
                  f"sequences file, e.g. {missing[:3]}.")
        sequences = {m: sequences[m] for m in wanted if m in sequences}
        if not sequences:
            raise ValueError(
                "None of the requested models were found in the sequences file."
            )

    organism = block.variables.get(organismVariable.id) or "Escherichia coli general"
    device = block.variables.get(deviceVariable.id) or "cpu"
    conda_env = block.variables.get(condaEnvVariable.id) or None
    add_start_methionine = block.variables.get(startMethionineVariable.id, True)

    print(f"Codon-optimising {len(sequences)} sequence(s) for '{organism}' on {device}...")
    if device == "cuda":
        print("Note: 'cuda' requires a GPU visible to the CodonTransformer environment.")

    try:
        # pylint: disable=import-outside-toplevel
        from bioprospecting.expression import codonTransformer
        # pylint: enable=import-outside-toplevel
    except ImportError as error:
        raise ValueError(
            "The 'bioprospecting' library is required by this block but is not "
            f"importable in the Horus backend ({error})."
        ) from error

    # bioprospecting activates the environment with 'conda activate', which is missing on
    # micromamba hosts. Put the environment's bin directory first on PATH instead and let
    # it run a plain 'python'.
    env_bin = _resolve_environment_bin(conda_env)
    python = _assert_codon_transformer(env_bin, conda_env)
    print(f"Running CodonTransformer with {python}")
    if env_bin:
        os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")

    codon_sequences = codonTransformer.computeCodonTransformerSequences(
        sequences,
        organism=organism,
        device=device,
        codon_transformer_conda_env=None,
        add_start_methionine=add_start_methionine,
    )

    codon_sequences = {str(k): str(v) for k, v in dict(codon_sequences).items()}

    missing = [m for m in sequences if m not in codon_sequences]
    if missing:
        print(f"Warning: no DNA returned for {len(missing)} model(s), "
              f"e.g. {missing[:3]}.")

    dna_json_output = "codon_optimised_dna.json"
    with open(dna_json_output, "w") as jf:
        json.dump(codon_sequences, jf, indent=2)

    dna_fasta_output = "codon_optimised_dna.fasta"
    with open(dna_fasta_output, "w") as of:
        for model, dna in codon_sequences.items():
            of.write(f">{model}\n{dna}\n")

    table_output = "codon_optimised_sequences.csv"
    with open(table_output, "w", newline="") as cf:
        writer = csv.DictWriter(
            cf,
            fieldnames=["model", "protein_length", "dna_length", "protein", "dna"],
        )
        writer.writeheader()
        for model, dna in codon_sequences.items():
            protein = sequences.get(model, "")
            writer.writerow(
                {
                    "model": model,
                    "protein_length": len(protein),
                    "dna_length": len(dna),
                    "protein": protein,
                    "dna": dna,
                }
            )

    print(f"Optimised {len(codon_sequences)} sequence(s).")

    rows = "".join(
        f"<tr><td>{model}</td><td>{len(sequences.get(model, ''))}</td>"
        f"<td>{len(dna)}</td>"
        f"<td style='font-family:monospace;font-size:11px;text-align:left'>"
        f"{dna[:60]}{'&hellip;' if len(dna) > 60 else ''}</td></tr>"
        for model, dna in codon_sequences.items()
    )
    html = f"""
    <style>
      table.eapm-table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
      table.eapm-table th, table.eapm-table td {{
          border: 1px solid #ddd; padding: 4px 8px; text-align: right;
      }}
      table.eapm-table th {{ background: #f2f2f2; }}
    </style>
    <div style="font-family:sans-serif">
      <h3>Codon-optimised sequences</h3>
      <p>{len(codon_sequences)} sequences &middot; host <b>{organism}</b></p>
      <div style="max-height:460px; overflow:auto">
        <table class="eapm-table">
          <tr><th>Model</th><th>Protein (aa)</th><th>DNA (nt)</th><th>DNA (start)</th></tr>
          {rows}
        </table>
      </div>
    </div>
    """
    Extensions().loadHTML(html, title="Codon-optimised sequences")

    block.setOutput(dnaFastaFile.id, dna_fasta_output)
    block.setOutput(dnaJsonFile.id, dna_json_output)
    block.setOutput(tableFile.id, table_output)


codonTransformerBlock = PluginBlock(
    category="Design & Prediction",
    name="CodonTransformer",
    id="codon_transformer",
    description="Codon-optimise protein sequences for expression in a host organism "
    "using CodonTransformer.",
    inputs=[sequencesFile, modelsFile],
    variables=[
        organismVariable,
        deviceVariable,
        condaEnvVariable,
        startMethionineVariable,
    ],
    outputs=[dnaFastaFile, dnaJsonFile, tableFile],
    action=compute_codon_sequences,
)
