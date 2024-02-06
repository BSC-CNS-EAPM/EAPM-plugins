"""
Module containing the Mafft block for the EAPM plugin
"""

import os
import subprocess
from HorusAPI import PluginVariable, PluginBlock, VariableTypes


# ==========================#
# Variable inputs
# ==========================#
sequences_mafft = PluginVariable(
    name="Sequences",
    id="sequences",
    description="The sequences for the msa",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# ==========================#
# Variable outputs
# ==========================#
output_mafft = PluginVariable(
    name="Output",
    id="output",
    description="The fasta output for the msa",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

def initialMafft(block: PluginBlock):
    sequences = block.inputs.get("sequences", None)
    output = block.outputs.get("output", None)
    
    try:
        subprocess.run(["mafft", "--version"], check=True)
    except Exception as exc:
        raise Exception(
            "MAFFT is not installed in the selected machine. Please install MAFFT in order to align the PDBs"
        ) 
    
    subprocess.run(["mafft", "--thread", "-1", "--auto", sequences, ">", output])
    
    block.setOutput("output", output)




