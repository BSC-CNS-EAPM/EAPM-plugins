"""
Module containing the AlphaFold block for the EAPM plugin
"""

import datetime
import os
import shutil
import subprocess
import tarfile
from HorusAPI import PluginVariable, PluginBlock, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
fasta_fileAF = PluginVariable(
    name="Fasta file",
    id="fasta_file",
    description="The input fasta file.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# ==========================#
# Variables
# ==========================#
outputAF = PluginVariable(
    name="Alphafold simulation folder",
    id="folder_name",
    description="The name of the folder where the simulation will be stored.",
    type=VariableTypes.STRING,
    defaultValue="alphafold",
)

# Output variables
outputJobsVariable = PluginVariable(
    id="output_jobs",
    name="BSC jobs",
    description="Output jobs from the docking calculation",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_jobs"],
)


# Alphafold action block
def initialAlphafold(block: PluginBlock):
    """
    Initial action of the block. It prepares the simulation and sends it to the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    # Loading plugin variables
    fastaFile = block.inputs.get("fasta_file", "None")
    if fastaFile == "None":
        raise Exception("No fasta file provided.")

    folderName = block.variables.get("folder_name", "alphafold")

    output_models_folder = os.path.join(os.getcwd(), folderName, "output_models")

    import prepare_proteins

    print("Loading fasta files...")

    sequences = prepare_proteins.sequenceModels(fastaFile)

    print("Setting up AlphaFold...")

    jobs = sequences.setUpAlphaFold(folderName)

    output_jobs = {
        "program": "alphafold",
        "jobs": jobs,
        "results_data": {
            "output_models": output_models_folder,
        },
    }

    print("Jobs created")

    block.setOutput("output_jobs", output_jobs)


alphafoldBlock = PluginBlock(
    name="Alphafold",
    description="Run Alphafold. (For cte_power, marenostrum and minotauro clusters or local)",
    action=initialAlphafold,
    variables=[],
    inputs=[fasta_fileAF],
    outputs=[outputJobsVariable],
)
