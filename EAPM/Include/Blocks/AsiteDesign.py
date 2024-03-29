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
    name="Input yaml",
    id="input_yaml",
    description=" Path to the input file yaml.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["yaml", "yml"],
)
inputPDBAsite = PluginVariable(
    name="Input PDB",
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
    id="folder_name",
    description="The name of the folder where the simulation will be stored.",
    type=VariableTypes.STRING,
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
    defaultValue=1,
)

queue = PluginVariable(
    name="Cluster queue",
    id="partition",
    description="The queue for the simulation",
    type=VariableTypes.STRING,
    defaultValue="bsc_ls",
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
    input_params = block.inputs.get("input_params", None)
    input_pdb = block.inputs.get("input_pdb", None)
    cpus = block.variables.get("cpus", 0)
    container = block.variables.get("container", None)
    output_file = input_yaml.rstrip(".yaml").split("/")[-1] + ".out"
    
    # copiar pdb and params to output folder
    subprocess.run(["cp", input_yaml, os.getcwd()], check=True)
    subprocess.run(["cp", input_pdb,  os.getcwd()], check=True)
    subprocess.run(["cp", "-r", input_params,  os.getcwd()], check=True)
    
    input_yaml = input_yaml.split("/")[-1]
    
    cluster = "local"
    if block.remote.name != "local":
        cluster = block.remote.host
    
    if "mn" in cluster:
        job = f"mpirun -np {cpus} python -m ActiveSiteDesign {input_yaml} > {output_file}"
    elif cluster == "local":
        if container is None:
            job = f"mpirun -np {cpus} python -m ActiveSiteDesign {input_yaml} > {output_file}"
        else:
            if cpus == 0:
                job = f"singularity exec {container} python -m ActiveSiteDesign {input_yaml} > {output_file}"
            else:
                job = f"mpirun -np {cpus} singularity exec {container} python -m ActiveSiteDesign {input_yaml} > {output_file}"    
    else:
        raise Exception("AsiteDesign can only be run on Marenostrum or local")
    
    
    from utils import launchCalculationAction

    launchCalculationAction(block, [job], "asitedesign", modulePurge=True)
    
           
def finalAsiteAction(block: SlurmBlock):
    from utils import downloadResultsAction
    downloaded_path = downloadResultsAction(block)

    resultsFolder = block.outputs["folder_name"]

    output_folder = os.path.join(downloaded_path, resultsFolder)

    block.setOutput(outputFolderAsite.id, output_folder)   


from utils import BSC_JOB_VARIABLES

asiteDesignBlock = SlurmBlock(
    name="AsiteDesign",
    description="Run AsiteDesign. (For local or marenostrum)",
    initialAction=initialAsite,
    finalAction=finalAsiteAction,
    variables=BSC_JOB_VARIABLES + [cpusAsite, containerAsite],
    inputs=[inputYamlAsite, inputPDBAsite, inputParamsAsite],
    outputs=[outputFolderAsite],
)

