"""
Module containing the rbcavity block for the EAPM plugin
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
outputLog = PluginVariable(
    name="Output log",
    id="output_log",
    description="The output log file.",
    type=VariableTypes.FILE,
    defaultValue="parameter_file.log",
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

    if block.remote.name != "Local":
        raise Exception("This block is only available for local execution.")

    # Loading plugin variables
    input_PRMfile = block.inputs.get(inputPRMFile.id, None)
    output_log = block.outputs.get(outputLog.id, "parameter_file.log")

    # rbcavity -was -d -r parameter_file.prm > parameter_file.log
    command = "rbcavity "
    if block.variables.get("was", True):
        command += "-was "
    if block.variables.get("dump_insight", False):
        command += "-d "
    command += f"-r {input_PRMfile} > {output_log}"

    print("Setting output of block to the results directory...")

    # subprocess the command
    import subprocess

    completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Get the output and error
    output = completed_process.stdout
    error = completed_process.stderr

    # Set the output
    block.setOutput(outputLog.id, output_log)


rbCavityBlock = PluginBlock(
    name="Rbcavity",
    description="Calculate docking cavities. (For local)",
    action=initialRbcavity,
    variables=[
        was,
        dumpInsight,
    ],
    inputs=[inputPRMFile],
    outputs=[outputLog],
)
