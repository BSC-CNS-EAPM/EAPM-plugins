"""
Module containing the rbdock block for the EAPM plugin
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
inputLigand = PluginVariable(
    name="Ligand SD file",
    id="input_ligand",
    description="The input ligand SD file.",
    type=VariableTypes.FILE,
    allowedValues=["sd"],
)

# ==========================#
# Variable outputs
# ==========================#
outputFile = PluginVariable(
    name="Output File",
    id="output_file",
    description="The output file with the ligand docked.",
    type=VariableTypes.FILE,
    defaultValue="parameter_file",
)


##############################
#       Other variables      #
##############################
protoPrmFile = PluginVariable(
    name="proto Prm File",
    id="proto_prm_file",
    description="The docking protocol parameter file.",
    type=VariableTypes.FILE,
    defaultValue="dock.prm",
)
nRuns = PluginVariable(
    name="nRuns",
    id="n_runs",
    description="Number of runs/ligand (default=1).",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)
allH = PluginVariable(
    name="allH",
    id="all_h",
    description="Keep all hydrogens, read all hydrogens present (default=polar hydrogens only).",
    type=VariableTypes.BOOLEAN,
    defaultValue=False,
)


# Align action block
def initialRbdock(block: PluginBlock):

    import os

    if block.remote.name != "Local":
        raise Exception("This block is only available for local execution.")

    # Loading plugin variables
    input_PRMfile = block.inputs.get(inputPRMFile.id, None)
    if input_PRMfile is None:
        raise Exception("No parameter file provided.")
    if not os.path.exists(input_PRMfile):
        raise Exception("Parameter file does not exist.")
    input_ligand = block.inputs.get(inputLigand.id, None)
    out = "output_dock"
    if input_ligand is None:
        raise Exception("No ligand file provided.")
    if not os.path.exists(input_ligand):
        raise Exception("Ligand file does not exist.")
    else:
        out = os.path.basename(input_ligand).split(".")[0] + "_out"
    output_file = block.outputs.get(outputFile.id, out)

    # rbcavity -was -d -r parameter_file.prm > parameter_file.log
    command = f"rbdock -i {input_ligand} -o {output_file} -r {input_PRMfile} "
    if block.variables.get("proto_prm_file", None) is not None:
        command += f"-p {block.variables.get('proto_prm_file')} "
    if block.variables.get("n_runs", None) is not None:
        command += f"-n {block.variables.get('n_runs')} "
    if block.variables.get("all_h", False):
        command += "-allH "

    print("Setting output of block to the results directory...")

    # Subprocess the command
    import subprocess

    completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Get the output and error
    output = completed_process.stdout
    # Save the output and error
    with open(f"{output_file}.out", "w") as f:
        f.write(output)
    error = completed_process.stderr
    with open(f"{output_file}.err", "w") as f:
        f.write(error)

    # Set the output
    block.setOutput(outputFile.id, output_file)


rbDockBlock = PluginBlock(
    name="Rbdock",
    description="Calculate the docking. (For local)",
    action=initialRbdock,
    variables=[
        protoPrmFile,
        nRuns,
        allH,
    ],
    inputs=[inputPRMFile, inputLigand],
    outputs=[outputFile],
)
