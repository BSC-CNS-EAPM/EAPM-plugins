"""
Module containing the sequence-length distribution block for the EAPM plugin.

Reproduces the UPO notebook's length check (cell 26): a histogram (with optional
KDE) of the selected sequences' lengths, annotated with the mean length.

Besides the plot, the block writes the per-model residue counts. The notebook
recomputes those repeatedly as::

    residue_lengths = {model: len(seq) for model, seq in models.sequences.items()}

to normalise Rosetta scores per residue, so emitting them once here lets the
downstream analysis blocks consume them directly.
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
sequencesFile = PluginVariable(
    id="sequences_file",
    name="Sequences file",
    description="Sequences to measure. A FASTA file or a JSON {name: sequence} mapping.",
    type=VariableTypes.FILE,
    allowedValues=["fasta", "fa", "faa", "json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
binsVariable = PluginVariable(
    id="bins",
    name="Bins",
    description="Number of histogram bins.",
    type=VariableTypes.INTEGER,
    defaultValue=30,
)
kdeVariable = PluginVariable(
    id="kde",
    name="Show KDE",
    description="Overlay a kernel density estimate on the histogram.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)
titleVariable = PluginVariable(
    id="title",
    name="Plot title",
    description="Title shown on the plot.",
    type=VariableTypes.STRING,
    defaultValue="Protein Sequence Length Distribution",
)

# ==========================#
# Variable outputs
# ==========================#
plotFile = PluginVariable(
    id="plot_file",
    name="Length histogram",
    description="PNG histogram of the sequence length distribution.",
    type=VariableTypes.FILE,
    allowedValues=["png"],
)
lengthsFile = PluginVariable(
    id="lengths_file",
    name="Sequence lengths",
    description="JSON mapping each model to its sequence length (residue count).",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
statsFile = PluginVariable(
    id="stats_file",
    name="Length statistics",
    description="JSON with count, mean, median, standard deviation, min and max length.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)


def _read_sequences(sequences_path: str) -> dict:
    """Read sequences from a FASTA or JSON file into a ``{name: sequence}`` dict."""
    import json

    if sequences_path.lower().endswith(".json"):
        with open(sequences_path) as jf:
            data = json.load(jf)
        if not isinstance(data, dict):
            raise ValueError(
                "The JSON sequences file must map sequence names to sequences "
                "(a plain list of names carries no length information)."
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


def plot_length_distribution(block: PluginBlock):
    """
    Plot the sequence length distribution and export per-model lengths.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import base64
    import json
    import statistics

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    # pylint: enable=import-outside-toplevel

    sequences_path = block.inputs.get(sequencesFile.id, None)
    if not sequences_path:
        raise ValueError("A valid sequences file (FASTA or JSON) must be provided.")

    bins = int(block.variables.get(binsVariable.id, 30) or 30)
    show_kde = block.variables.get(kdeVariable.id, True)
    title = block.variables.get(titleVariable.id) or "Protein Sequence Length Distribution"

    sequences = _read_sequences(sequences_path)

    lengths = {name: len(seq) for name, seq in sequences.items()}
    empty = [name for name, length in lengths.items() if length == 0]
    if empty:
        print(f"Warning: {len(empty)} sequence(s) are empty and will skew the "
              f"distribution, e.g. {empty[:3]}.")

    values = list(lengths.values())
    if not values:
        raise ValueError("No sequences to measure.")

    mean_length = sum(values) / len(values)
    stats = {
        "count": len(values),
        "mean": mean_length,
        "median": statistics.median(values),
        "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
        "min": min(values),
        "max": max(values),
    }

    print(f"{stats['count']} sequences: mean {stats['mean']:.1f}, "
          f"median {stats['median']}, range {stats['min']}-{stats['max']} residues.")

    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))

    # A KDE needs some spread; fall back to a plain histogram when it cannot be fit.
    try:
        sns.histplot(values, bins=bins, kde=bool(show_kde),
                     color="skyblue", edgecolor="black")
    except Exception as error:  # pylint: disable=broad-except
        print(f"Could not draw the KDE ({error}); falling back to a plain histogram.")
        plt.clf()
        sns.histplot(values, bins=bins, kde=False, color="skyblue", edgecolor="black")

    plt.axvline(mean_length, color="red", linestyle="--", linewidth=2,
                label=f"Mean = {mean_length:.1f}")
    plt.title(title, fontsize=16)
    plt.xlabel("Sequence Length (amino acids)", fontsize=12)
    plt.ylabel("Count", fontsize=12)
    plt.legend()
    plt.tight_layout()

    plot_output = "sequence_length_distribution.png"
    plt.savefig(plot_output, dpi=150)
    plt.close()

    with open(plot_output, "rb") as img:
        encoded = base64.b64encode(img.read()).decode("ascii")
    html = (
        f'<div style="text-align:center">'
        f"<h3>{title}</h3>"
        f'<img src="data:image/png;base64,{encoded}" style="max-width:100%"/>'
        f"<p>{stats['count']} sequences &middot; mean {stats['mean']:.1f} &middot; "
        f"range {stats['min']}-{stats['max']} residues</p>"
        f"</div>"
    )
    Extensions().loadHTML(html, title="Sequence length distribution")

    lengths_output = "sequence_lengths.json"
    with open(lengths_output, "w") as jf:
        json.dump(lengths, jf, indent=2)

    stats_output = "sequence_length_stats.json"
    with open(stats_output, "w") as jf:
        json.dump(stats, jf, indent=2)

    block.setOutput(plotFile.id, plot_output)
    block.setOutput(lengthsFile.id, lengths_output)
    block.setOutput(statsFile.id, stats_output)


sequenceLengthDistributionBlock = PluginBlock(
    category="Clustering & Selection",
    name="Sequence Length Distribution",
    id="sequence_length_distribution",
    description="Plot the length distribution of a set of sequences and export the "
    "per-model residue counts.",
    inputs=[sequencesFile],
    variables=[binsVariable, kdeVariable, titleVariable],
    outputs=[plotFile, lengthsFile, statsFile],
    action=plot_length_distribution,
)
