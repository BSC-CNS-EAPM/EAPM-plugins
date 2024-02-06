"""
Module containing the Asitedesign block for the EAPM plugin
"""

import subprocess
import os
from HorusAPI import PluginVariable, VariableTypes, SlurmBlock, VariableList

# ==========================#
# Variable inputs
# ==========================#
inputYamlAsite= PluginVariable(
    name="PDB reference",
    id="input_yaml",
    description=" Path to the input file yaml.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["yaml", "yml"],
)
inputPDBAsite = PluginVariable(
    name="PDB reference",
    id="input_pdb",
    description=" Path to the input file PDB.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["pdb"],
)
inputParamsAsite = PluginVariable(
    name="Parameters",
    id="input_params",
    description=" Path to the input file parameters folder.",
    type=VariableTypes.FOLDER,
    defaultValue=None,
)

# ==========================#
# Variable outputs
# ==========================#
outputFolderAsite = PluginVariable(
    name="Asite simulation folder",
    id="folder_output",
    description="The name of the folder where the simulation will be stored.",
    type=VariableTypes.Folder,
    defaultValue="AsiteDesign",
)

##############################
#       Other variables      #
##############################
cpusAsite = PluginVariable(
    name="CPUs",
    id="cpus",
    description="Number of CPUs to use for the calculation.",
    type=VariableTypes.INTEGER,
    defaultValue=0,
)
folder_nameAsite = PluginVariable(
    name="Asite simulation folder",
    id="folder_name",
    description="The name of the folder where the simulation will be stored.",
    type=VariableTypes.STRING,
    defaultValue="Asite",
)
containerAsite = PluginVariable(
    name="Container",
    id="container",
    description="If you are launching the block in a container. The container to use.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
  
def initialAsite(block: SlurmBlock):
    # Get the input variables
    input_yaml = block.inputs.get("input_yaml", None)
    cpus = block.variables.get("cpus", 0)
    output = block.outputs.get("output_Asite", None)
    container = block.variables.get("container", None)
    output_file = input_yaml.rstrip(".yaml").split("/")[-1] + ".out"
    
    if "mn" in block.remote.host:
        job = f"mpirun -np {cpus} python -m ActiveSiteDesign {input_yaml} > {output_file}"
        
    else:
        if container is None:
            job = f"mpirun -np {cpus} python -m ActiveSiteDesign {input_yaml} > {output_file}"
        else:
            if cpus == 0:
                job = f"singularity exec {container} python -m ActiveSiteDesign {input_yaml} > {output_file}"
            else:
                job = f"mpirun -np {cpus} singularity exec {container} python -m ActiveSiteDesign {input_yaml} > {output_file}"    
    
    from utils import launchCalculationAction

    launchCalculationAction(block, job)
    
           
def finalAsiteAction(block: SlurmBlock):
    from utils import downloadResultsAction
    downloaded_path = downloadResultsAction(block)

    resultsFolder = block.output["folder_output"]

    output_folder = os.path.join(downloaded_path, resultsFolder)

    block.setOutput(outputFolderAsite.id, output_folder)   

asiteDesignBlock = SlurmBlock(
    name="AsiteDesign",
    description="Run AsiteDesign. (For local or marenostrum)",
    initialAction=initialAsite,
    finalAction=finalAsiteAction,
    variables=[cpusAsite, folder_nameAsite, containerAsite],
    inputs=[inputYamlAsite, inputPDBAsite, inputParamsAsite],
    outputs=[outputFolderAsite],
)

