"""
Module containing the HmmSearch block for the EAPM plugin as a local implementation
"""

import datetime
import os
import pyhmmer
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
sequenceDBVar = PluginVariable(
    id="sequence_db",
    name="Sequence DB",
    description="The sequence database to search",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/bsc72/ruite/enzyminer/data/reduced_data/uniref50.fasta",
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
exeVar = PluginVariable(
    id="hmmsearch_exe",
    name="Hmmsearch executable",
    description="The hmmsearch executable",
    type=VariableTypes.STRING,
    defaultValue="/gpfs/projects/bsc72/anarobles/HMM_proof/hmmer-3.3.2/src/hmmsearch", 
)
cpuVar = PluginVariable(
    id="hmmsearch_cpu",
    name="Hmmsearch CPU",
    description="The number of CPUs to use",
    type=VariableTypes.INTEGER,
    defaultValue=10,
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
    
    if "mn" not in block.remote.host:
        raise Exception("This block only works on Marenostrum.")
    
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
    
    exe = block.variables["hmmsearch_exe"]
    cpu = block.variables["hmmsearch_cpu"]
    evalue = block.variables["hmmsearch_evalue"]
    
    command = f"{exe} --cpu {cpu} -E {evalue} {input} {sequenceDBVar}"
    
    #! Finish the command execution in the remote (see if it's necessary to use a sbatch file)
    
    output = block.outputs.get("output", "output.domtbl")
    
        
    block.setOutput("outputVariable", output)
    
def finalAction():
    pass


hmmsearchBlock = SlurmBlock(
    name="HmmSearch",
    initialAction=runHmmSearch,
    finalAction=finalAction,
    description="Searches a sequence database with a given hmm",
    inputs=[hmmInput, sequenceDBVar],
    variables=[],
    outputs=[outputVariable],
)
