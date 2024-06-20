"""
Module containing the HmmSearch block for the EAPM plugin as a nord3 implementation
"""

import datetime
import os
from HorusAPI import SlurmBlock, PluginVariable, VariableTypes, Extensions


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

# ==========================#
# Variable
# ==========================#
sequenceDBVar = PluginVariable(
    id="sequence_db",
    name="Sequence DB",
    description="The sequence database to search",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/shared/public/AlphaFold/uniref90/uniref90.fa",
)
evalueVar = PluginVariable(
    id="hmmsearch_evalue",
    name="Hmmsearch evalue",
    description="The evalue to use",
    type=VariableTypes.FLOAT,
    defaultValue=0.001,
)


def runHmmSearch(block: SlurmBlock):

    input = block.inputs.get("input_hmm", None)

    if "nord3" not in block.remote.host:
        raise Exception("This block only works on Nord3.")

    if input is None:
        raise Exception("No input hmm provided")
    if not os.path.exists(input):
        raise Exception(f"The input hmm file does not exist: {input}")

    savedID_and_date = block.flow.savedID + "_" + str(datetime.datetime.now().timestamp())
    simRemoteDir = os.path.join(block.remote.workDir, savedID_and_date)
    block.extraData["remoteDir"] = simRemoteDir
    block.remote.remoteCommand(f"mkdir -p -v {simRemoteDir}")

    print(f"Created simulation folder in the remote at {simRemoteDir}")
    print("Sending data to the remote...")

    exe = block.config.get("hmmer_path", "hmmer")
    cpu = block.variables["hmmsearch_cpu"]
    evalue = block.variables["hmmsearch_evalue"]

    command = f"{exe} --cpu {cpu} -E {evalue} {input} {sequenceDBVar}"

    #! Finish the command execution in the remote (see if it's necessary to use a sbatch file)

    output = block.outputs.get("output", "output.domtbl")

    block.setOutput("outputVariable", output)


def finalAction():
    pass


from utils import BSC_JOB_VARIABLES

hmmsearchBlock = SlurmBlock(
    name="HmmSearch",
    initialAction=runHmmSearch,
    finalAction=finalAction,
    description="Searches a sequence database with a given hmm",
    inputs=[hmmInput],
    variables=BSC_JOB_VARIABLES + [sequenceDBVar, evalueVar],
    outputs=[outputVariable],
)
