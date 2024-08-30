"""
Module containing the rDock cavity block for the EAPM plugin
"""

from HorusAPI import PluginBlock, PluginVariable, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
inputPRMFile = PluginVariable(
    name="Parameter file",
    id="input_prm_file",
    description="The input '.prm' file.",
    type=VariableTypes.FILE,
    defaultValue="parameter_file.prm",
    allowedValues=["prm"],
)

# ==========================#
# Variable outputs
# ==========================#
cavityFile = PluginVariable(
    name="Cavity file",
    id="cavity_file",
    description="File required to perform an rDock docking.",
    type=VariableTypes.FILE,
    defaultValue=None,
)

##############################
#       Other variables      #
##############################
was = PluginVariable(
    name="Was",
    id="was",
    description="Write docking cavities (plus distance grid) to .as file.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
)
dumpInsight = PluginVariable(
    name="Dump Insight",
    id="dump_insight",
    description="Dump InsightII/PyMOL grids for each cavity for visualisation.",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)


# Align action block
def initialRbcavity(block: PluginBlock):

    import os

    if block.remote.name != "Local":
        raise Exception("This block is only available for local execution.")

    # Loading plugin variables
    input_PRMfile = block.inputs.get(inputPRMFile.id, None)
    output_log = "cavity_generation.log"

    path_output = os.path.dirname(input_PRMfile)

    # rbcavity -was -d -r parameter_file.prm > parameter_file.log
    command = "rbcavity "
    if block.variables.get("was", True):
        command += "-was "
    if block.variables.get("dump_insight", False):
        command += "-d "
    command += f"-r {input_PRMfile} > {os.path.join(path_output,output_log)}"

    # subprocess the command
    import subprocess

    completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Get the output and error
    output = completed_process.stdout
    error = completed_process.stderr

    # Save the output and error
    with open(f"{os.path.join(path_output,'cavity_grid.out')}", "w") as f:
        f.write(output)
    with open(f"{os.path.join(path_output,'cavity_grid.err')}", "w") as f:
        f.write(error)

    parameter_file_name = input_PRMfile.split(".")[0]

    # Set the output
    block.setOutput(cavityFile.id, os.path.join(path_output, f"{parameter_file_name}.as"))


rbCavityBlock = PluginBlock(
    name="rDockCavity",
    id="rbcavity",
    description="Calculate docking cavity for rDock. (For local)",
    action=initialRbcavity,
    variables=[
        was,
        dumpInsight,
    ],
    inputs=[inputPRMFile],
    outputs=[cavityFile],
)
