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


def convertMSA2HMM(block: PluginBlock):
    """
    Convert MSA to HMM
    """
    import os

    import pyhmmer

    # Loading plugin variables
    inputMSA = block.inputs.get("input_file_msa")
    if inputMSA is None:
        raise Exception("No input MSA provided")

    if not os.path.exists(inputMSA):
        raise Exception(f"The input MSA file does not exist: {inputMSA}")

    alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.MSAFile(inputMSA, digital=True, alphabet=alphabet) as msa_file:
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
    action=convertMSA2HMM,
)
