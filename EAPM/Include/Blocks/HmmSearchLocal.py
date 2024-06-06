"""
Module containing the HmmSearch block for the EAPM plugin as a local implementation
"""

from HorusAPI import Extensions, PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
hmmInput = PluginVariable(
    id="input_hmm",
    name="Hmm input",
    description="The input hmm",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["hmm"],
)
sequenceDBVar = PluginVariable(
    id="sequence_db",
    name="Sequence DB",
    description="The sequence database to search",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta", "faa"],
)

# ==========================#
# Variable outputs
# ==========================#
outputVariable = PluginVariable(
    id="output",
    name="Output File",
    description="Output of the Hmmsearch block",
    type=VariableTypes.FILE,
    allowedValues=["domtbl"],
)


def runHmmSearch(block: PluginBlock):
    import os

    import pyhmmer

    input = block.inputs.get("input_hmm", None)

    if input is None:
        raise Exception("No input hmm provided")

    if not os.path.exists(input):
        raise Exception(f"The input hmm file does not exist: {input}")

    try:
        with pyhmmer.plan7.HMMFile(input) as hmm_file:
            hmm = hmm_file.read()
    except Exception as e:
        raise Exception(f"Error reading the input hmm file: {e}")

    alphabet = pyhmmer.plan7.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)

    sequenceDB = block.inputs.get("sequence_db", None)

    if sequenceDB is None:
        raise Exception("No sequence database provided")

    if not os.path.exists(sequenceDB):
        raise Exception(f"The sequence database file does not exist: {sequenceDB}")

    try:
        with pyhmmer.easel.SequenceFile(sequenceDB, digital=True, alphabet=alphabet) as seq_file:
            hits = pipeline.search_hmm(hmm, seq_file)
    except Exception as e:
        raise Exception(f"Error searching the sequence database: {e}")

    output = block.outputs.get("output", "output.domtbl")

    with open(output, "wb") as f:
        hits.write(f, format="domains")

    block.setOutput("outputVariable", output)


hmmsearchLocalBlock = PluginBlock(
    name="HmmSearch Local",
    id="hmmsearch_local",
    description="Searches a sequence database with a given hmm",
    inputs=[hmmInput, sequenceDBVar],
    variables=[],
    outputs=[outputVariable],
    action=runHmmSearch,
)
