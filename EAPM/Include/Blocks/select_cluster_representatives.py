"""
Module containing the cluster-representative selection block for the EAPM plugin.

Given an MMseqs2 clustering (``{representative: [members]}``) and a per-model
score (e.g. ProteinMPNN mean scores), this block reproduces the selection logic
of the UPO notebook: rank clusters by the best score among their members, keep
the top-N clusters, and emit the best-scoring member of each as the selected
representative set.

By default lower scores are considered better (as with ProteinMPNN negative
log-likelihoods). Set 'Maximize score' to invert this.
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
clustersFile = PluginVariable(
    id="clusters_file",
    name="Clusters file",
    description="JSON file mapping each cluster representative to its members "
    "(output of the MMseqs2 Clustering block).",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
scoresFile = PluginVariable(
    id="scores_file",
    name="Scores file",
    description="Per-model scores. Either a JSON file mapping model name to score, "
    "or a CSV with a name column and a score column (e.g. 'mean_score').",
    type=VariableTypes.FILE,
    allowedValues=["json", "csv"],
)
sequencesFile = PluginVariable(
    id="sequences_file",
    name="Sequences file (optional)",
    description="Optional FASTA or JSON ({name: sequence}) file. If given, the "
    "sequences of the selected representatives are written to a FASTA output.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta", "fa", "faa", "json"],
)

# ==========================#
# Variables (parameters)
# ==========================#
topNVariable = PluginVariable(
    id="top_n_clusters",
    name="Top N clusters",
    description="Number of best-scoring clusters to keep (one representative each).",
    type=VariableTypes.INTEGER,
    defaultValue=100,
)
maximizeVariable = PluginVariable(
    id="maximize",
    name="Maximize score",
    description="If enabled, higher scores are better. By default (disabled) lower "
    "scores are better, as with ProteinMPNN scores.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)

# ==========================#
# Variable outputs
# ==========================#
selectedModelsFile = PluginVariable(
    id="selected_models_file",
    name="Selected models",
    description="JSON list with the selected representative model of each top cluster.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
selectedScoresFile = PluginVariable(
    id="selected_scores_file",
    name="Selected scores",
    description="JSON mapping each selected model to its score.",
    type=VariableTypes.FILE,
    allowedValues=["json"],
)
selectedSequencesFile = PluginVariable(
    id="selected_sequences_file",
    name="Selected sequences",
    description="FASTA with the sequences of the selected models (only if a "
    "sequences file was provided).",
    type=VariableTypes.FILE,
    allowedValues=["fasta"],
)
selectionTableFile = PluginVariable(
    id="selection_table_file",
    name="Selection table",
    description="CSV summarising each selected cluster: representative, best model, "
    "best score, cluster size and number of scored members.",
    type=VariableTypes.FILE,
    allowedValues=["csv"],
)


def _load_scores(scores_path: str) -> dict:
    """Read per-model scores from a JSON ({name: score}) or CSV file."""
    import csv
    import json

    if scores_path.lower().endswith(".json"):
        with open(scores_path) as jf:
            data = json.load(jf)
        if not isinstance(data, dict):
            raise ValueError("The JSON scores file must map model names to scores.")
        return {str(name): float(value) for name, value in data.items()}

    # CSV: detect the name column (first, or a labelled one) and the score column.
    with open(scores_path, newline="") as cf:
        rows = list(csv.reader(cf))
    if not rows:
        raise ValueError("The scores CSV file is empty.")

    header = [h.strip() for h in rows[0]]
    lower = [h.lower() for h in header]

    name_idx = 0
    for label in ("model", "description", "name", "id"):
        if label in lower:
            name_idx = lower.index(label)
            break

    score_idx = None
    for label in ("mean_score", "score", "avg_score"):
        if label in lower:
            score_idx = lower.index(label)
            break
    if score_idx is None:
        score_idx = 1 if len(header) > 1 else None
    if score_idx is None:
        raise ValueError("Could not find a score column in the scores CSV file.")

    scores: dict = {}
    for row in rows[1:]:
        if len(row) <= max(name_idx, score_idx):
            continue
        name = row[name_idx].strip()
        raw = row[score_idx].strip()
        if not name or not raw:
            continue
        try:
            scores[name] = float(raw)
        except ValueError:
            continue
    if not scores:
        raise ValueError("No valid (name, score) pairs were read from the scores CSV file.")
    return scores


def _load_sequences(sequences_path: str) -> dict:
    """Read sequences from a FASTA or JSON file into a ``{name: sequence}`` dict."""
    import json

    if sequences_path.lower().endswith(".json"):
        with open(sequences_path) as jf:
            data = json.load(jf)
        if not isinstance(data, dict):
            raise ValueError("The JSON sequences file must map sequence names to sequences.")
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
    return sequences


def select_representatives(block: PluginBlock):
    """
    Rank clusters by their best member score and select the top-N representatives.

    Args:
        block (PluginBlock): The PluginBlock object representing the current block.
    """
    # pylint: disable=import-outside-toplevel
    import csv
    import json
    import os

    # pylint: enable=import-outside-toplevel

    clusters_path = block.inputs.get(clustersFile.id, None)
    if not clusters_path or not os.path.isfile(clusters_path):
        raise ValueError("A valid clusters JSON file must be provided.")

    scores_path = block.inputs.get(scoresFile.id, None)
    if not scores_path or not os.path.isfile(scores_path):
        raise ValueError("A valid scores file (JSON or CSV) must be provided.")

    sequences_path = block.inputs.get(sequencesFile.id, None)

    top_n = block.variables.get(topNVariable.id, 100)
    maximize = block.variables.get(maximizeVariable.id, False)

    with open(clusters_path) as jf:
        clusters = json.load(jf)
    if not isinstance(clusters, dict):
        raise ValueError("The clusters file must map representatives to member lists.")

    scores = _load_scores(scores_path)

    inf = float("inf")

    def cluster_best_score(members):
        vals = [scores[m] for m in members if m in scores]
        if not vals:
            return None
        return max(vals) if maximize else min(vals)

    def sort_key(cluster_id):
        best = cluster_best_score(clusters[cluster_id])
        if best is None:
            return inf  # clusters with no scored members sink to the bottom
        return -best if maximize else best

    sorted_clusters = sorted(clusters.keys(), key=sort_key)
    top_clusters = sorted_clusters[: int(top_n)]

    selected_models: list = []
    selected_scores: dict = {}
    table_rows: list = []

    for cluster_id in top_clusters:
        members = clusters[cluster_id]
        scored_members = [m for m in members if m in scores]
        if not scored_members:
            continue

        if maximize:
            best_model = max(scored_members, key=lambda m: scores[m])
        else:
            best_model = min(scored_members, key=lambda m: scores[m])

        selected_models.append(best_model)
        selected_scores[best_model] = scores[best_model]
        table_rows.append(
            {
                "cluster_representative": cluster_id,
                "best_model": best_model,
                "best_score": scores[best_model],
                "cluster_size": len(members),
                "scored_members": len(scored_members),
            }
        )

    print(
        f"Ranked {len(clusters)} clusters; kept top {len(top_clusters)}; "
        f"selected {len(selected_models)} representatives."
    )

    # Write outputs into the flow working directory
    selected_models_output = "selected_models.json"
    with open(selected_models_output, "w") as jf:
        json.dump(selected_models, jf, indent=2)
    block.setOutput(selectedModelsFile.id, selected_models_output)

    selected_scores_output = "selected_scores.json"
    with open(selected_scores_output, "w") as jf:
        json.dump(selected_scores, jf, indent=2)
    block.setOutput(selectedScoresFile.id, selected_scores_output)

    selection_table_output = "selection_table.csv"
    with open(selection_table_output, "w", newline="") as cf:
        writer = csv.DictWriter(
            cf,
            fieldnames=[
                "cluster_representative",
                "best_model",
                "best_score",
                "cluster_size",
                "scored_members",
            ],
        )
        writer.writeheader()
        writer.writerows(table_rows)
    block.setOutput(selectionTableFile.id, selection_table_output)

    if sequences_path and os.path.isfile(sequences_path):
        sequences = _load_sequences(sequences_path)
        missing = [m for m in selected_models if m not in sequences]
        if missing:
            print(f"Warning: {len(missing)} selected models are missing from the "
                  f"sequences file (e.g. {missing[:3]}).")
        selected_sequences_output = "selected_sequences.fasta"
        with open(selected_sequences_output, "w") as of:
            for model in selected_models:
                if model in sequences:
                    of.write(f">{model}\n{sequences[model]}\n")
        block.setOutput(selectedSequencesFile.id, selected_sequences_output)


selectClusterRepresentativesBlock = PluginBlock(
    category="Clustering & Selection",
    name="Select Cluster Representatives",
    id="select_cluster_representatives",
    description="Rank MMseqs2 clusters by their best member score and select the "
    "best-scoring representative of the top-N clusters.",
    inputs=[clustersFile, scoresFile, sequencesFile],
    variables=[topNVariable, maximizeVariable],
    outputs=[
        selectedModelsFile,
        selectedScoresFile,
        selectedSequencesFile,
        selectionTableFile,
    ],
    action=select_representatives,
)
