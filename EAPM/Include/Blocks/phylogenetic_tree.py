"""
Module containing the phylogenetic tree block for the EAPM plugin.

Reproduces the UPO notebook's tree construction (cells 75-78)::

    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(msa)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)          # or .upgma(dm)
    # terminals renamed to "seq_id (organism)"
    draw_ascii(tree)

The block takes an alignment (e.g. the MAFFT block's output), builds a
neighbour-joining or UPGMA tree from a substitution-matrix distance, and renders
it both as ASCII and as a figure. Terminal labels can be annotated with an
organism (or any other) name from a JSON mapping.

Pure Biopython, so there is no external tool dependency.
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
msaFile = PluginVariable(
    id="msa_file",
    name="Alignment file",
    description="Multiple sequence alignment (FASTA, Clustal, Phylip, Stockholm...). "
    "The format is guessed from the extension unless one is set below.",
    type=VariableTypes.FILE,
    allowedValues=["fasta", "fa", "aln", "clustal", "maf", "phy", "sto"],
)
labelsFile = PluginVariable(
    id="labels_file",
    name="Labels file (optional)",
    description="Optional JSON {sequence_id: label} used to annotate the leaves, "
    "e.g. an organism per model. Leaves become 'id (label)'.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
methodVariable = PluginVariable(
    id="method",
    name="Tree method",
    description="Tree construction method.",
    type=VariableTypes.STRING,
    defaultValue="nj",
    allowedValues=["nj", "upgma"],
)
modelVariable = PluginVariable(
    id="model",
    name="Distance model",
    description="Substitution model used for the pairwise distances.",
    type=VariableTypes.STRING,
    defaultValue="blosum62",
    allowedValues=[
        "identity",
        "blosum62",
        "blosum45",
        "blosum80",
        "blosum90",
        "pam250",
        "pam70",
        "pam30",
        "blastp",
    ],
)
formatVariable = PluginVariable(
    id="alignment_format",
    name="Alignment format",
    description="Force the alignment format instead of guessing from the extension.",
    type=VariableTypes.STRING,
    defaultValue="auto",
    allowedValues=["auto", "fasta", "clustal", "phylip", "stockholm", "maf"],
)

# ==========================#
# Variable outputs
# ==========================#
newickFile = PluginVariable(
    id="newick_file",
    name="Newick tree",
    description="The tree in Newick format.",
    type=VariableTypes.FILE,
    allowedValues=["nwk"],
)
asciiFile = PluginVariable(
    id="ascii_file",
    name="ASCII tree",
    description="Text rendering of the tree.",
    type=VariableTypes.FILE,
    allowedValues=["txt"],
)
plotFile = PluginVariable(
    id="plot_file",
    name="Tree plot",
    description="PNG rendering of the tree.",
    type=VariableTypes.FILE,
    allowedValues=["png"],
)

# Extension -> Biopython AlignIO format
_FORMAT_BY_EXTENSION = {
    ".fasta": "fasta",
    ".fa": "fasta",
    ".faa": "fasta",
    ".aln": "clustal",
    ".clustal": "clustal",
    ".maf": "clustal",
    ".phy": "phylip",
    ".phylip": "phylip",
    ".sto": "stockholm",
    ".stockholm": "stockholm",
}


def _read_alignment(path: str, forced_format: str):
    """Read an alignment, guessing the format from the extension when needed."""
    import os

    from Bio import AlignIO

    if forced_format and forced_format != "auto":
        return AlignIO.read(path, forced_format)

    extension = os.path.splitext(path)[1].lower()
    guess = _FORMAT_BY_EXTENSION.get(extension)

    candidates = [guess] if guess else []
    # Fall back to trying the common formats so an odd extension still works.
    candidates += [f for f in ("fasta", "clustal", "phylip", "stockholm") if f != guess]

    errors = []
    for candidate in candidates:
        try:
            alignment = AlignIO.read(path, candidate)
            print(f"Read the alignment as '{candidate}'.")
            return alignment
        except Exception as error:  # pylint: disable=broad-except
            errors.append(f"{candidate}: {error}")

    raise ValueError(
        "Could not read the alignment. Tried:\n  " + "\n  ".join(errors)
    )


def build_tree(block: PluginBlock):
    """
    Build a phylogenetic tree from an alignment and render it.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import base64
    import io
    import json
    import os

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

    # pylint: enable=import-outside-toplevel

    msa_path = block.inputs.get(msaFile.id, None)
    if not msa_path or not os.path.isfile(msa_path):
        raise ValueError("A valid alignment file must be provided.")

    method = (block.variables.get(methodVariable.id) or "nj").lower()
    if method not in ("nj", "upgma"):
        raise ValueError("Method must be 'nj' or 'upgma'.")

    model = block.variables.get(modelVariable.id) or "blosum62"
    alignment_format = block.variables.get(formatVariable.id) or "auto"

    alignment = _read_alignment(msa_path, alignment_format)
    if len(alignment) < 3:
        raise ValueError(
            f"A tree needs at least 3 sequences; the alignment has {len(alignment)}."
        )
    print(f"Alignment: {len(alignment)} sequences x {alignment.get_alignment_length()} columns.")

    labels = {}
    labels_path = block.inputs.get(labelsFile.id, None)
    if labels_path and os.path.isfile(labels_path):
        with open(labels_path) as jf:
            labels = {str(k): str(v) for k, v in json.load(jf).items()}
        print(f"Loaded {len(labels)} leaf labels.")

    print(f"Computing {model} distances...")
    calculator = DistanceCalculator(model)
    distance_matrix = calculator.get_distance(alignment)

    print(f"Building the {method.upper()} tree...")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix) if method == "nj" else constructor.upgma(
        distance_matrix
    )

    # Biopython names internal nodes "Inner1", "Inner2", ...; hide them in the render.
    for clade in tree.get_nonterminals():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = None

    if labels:
        annotated = 0
        for terminal in tree.get_terminals():
            label = labels.get(terminal.name)
            if label:
                terminal.name = f"{terminal.name} ({label})"
                annotated += 1
        print(f"Annotated {annotated}/{len(tree.get_terminals())} leaves.")

    newick_output = "phylogenetic_tree.nwk"
    Phylo.write(tree, newick_output, "newick")

    ascii_buffer = io.StringIO()
    Phylo.draw_ascii(tree, file=ascii_buffer)
    ascii_text = ascii_buffer.getvalue()

    ascii_output = "phylogenetic_tree.txt"
    with open(ascii_output, "w") as tf:
        tf.write(ascii_text)

    height = max(4.0, 0.28 * len(tree.get_terminals()))
    figure, axis = plt.subplots(figsize=(11, height))
    Phylo.draw(tree, axes=axis, do_show=False)
    axis.set_title(f"{method.upper()} tree ({model} distances)")
    figure.tight_layout()
    plot_output = "phylogenetic_tree.png"
    figure.savefig(plot_output, dpi=150)
    plt.close(figure)

    with open(plot_output, "rb") as img:
        encoded = base64.b64encode(img.read()).decode("ascii")

    escaped = ascii_text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    html = f"""
    <div style="font-family:sans-serif">
      <h3>Phylogenetic tree</h3>
      <p>{len(alignment)} sequences &middot; {method.upper()} &middot; {model} distances</p>
      <div style="overflow-x:auto; margin-bottom:24px">
        <img src="data:image/png;base64,{encoded}" style="max-width:100%"/>
      </div>
      <h3>ASCII</h3>
      <pre style="font-size:11px; line-height:1.2; overflow:auto; max-height:420px;
                  background:#f7f7f7; padding:12px">{escaped}</pre>
    </div>
    """
    Extensions().loadHTML(html, title="Phylogenetic tree")

    block.setOutput(newickFile.id, newick_output)
    block.setOutput(asciiFile.id, ascii_output)
    block.setOutput(plotFile.id, plot_output)


phylogeneticTreeBlock = PluginBlock(
    category="Sequence Alignment",
    name="Phylogenetic Tree",
    id="phylogenetic_tree",
    description="Build a neighbour-joining or UPGMA tree from a multiple sequence "
    "alignment and render it as Newick, ASCII and a figure.",
    inputs=[msaFile, labelsFile],
    variables=[methodVariable, modelVariable, formatVariable],
    outputs=[newickFile, asciiFile, plotFile],
    action=build_tree,
)
