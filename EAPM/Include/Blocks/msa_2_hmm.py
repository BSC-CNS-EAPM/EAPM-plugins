"""
Module containing the MSA2HMM block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableGroup, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
inputFileMSA = PluginVariable(
    id="input_file_msa",
    name="Input File MSA",
    description="The input file MSA",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)
inputMSAVar = PluginVariable(
    id="input_msa",
    name="Input MSA",
    description="The input MSA",
    type=VariableTypes.CUSTOM,
    defaultValue=None,
    allowedValues=["msa"],
)

# ==========================#
# Variable outputs
# ==========================#
hmmOutput = PluginVariable(
    id="output_hmm",
    name="Output HMM",
    description="Output file of the Hmmer block",
    type=VariableTypes.FILE,
    allowedValues=["hmm"],
)

# ==========================#
# Variables
# ==========================#


def convert_msa_2_hmm(block: PluginBlock):
    """
    Convert MSA to HMM
    """
    # pylint: disable=import-outside-toplevel
    import os

    import pyhmmer
    from pyhmmer.easel import MSAFile

    # pylint: enable=import-outside-toplevel
    # Loading plugin variables
    input_msa = block.inputs.get(inputFileMSA.id, None)
    if input_msa is None:
        raise ValueError("No input MSA provided")

    if not os.path.exists(input_msa):
        raise ValueError(f"The input MSA file does not exist: {input_msa}")

    alphabet = pyhmmer.easel.Alphabet.amino()

    with MSAFile(input_msa, digital=True, alphabet=alphabet) as msa_file:
        msa = msa_file.read()
        msa.name = b"input_msa"

    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)

    output = "output.hmm"
    with open(output, "wb") as output_file:
        hmm.write(output_file)

    block.setOutput("output_hmm", output)


convertMSAToHMMBlock = PluginBlock(
    name="MSA to HMM",
    id="msa_to_hmm",
    description="Convert MSA files to HMM",
    inputGroups=[
        VariableGroup(
            id=inputFileMSA.id,
            name=inputFileMSA.name,
            description=inputFileMSA.description,
            variables=[inputFileMSA],
        ),
        VariableGroup(
            id=inputMSAVar.id,
            name=inputMSAVar.name,
            description=inputMSAVar.description,
            variables=[inputMSAVar],
        ),
    ],
    variables=[],
    outputs=[hmmOutput],
    action=convert_msa_2_hmm,
)
