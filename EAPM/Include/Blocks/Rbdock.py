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
cavityFile = PluginVariable(
    name="Cavity File",
    id="cavity_file",
    description="Cavity file of the docking volume.",
    type=VariableTypes.FILE,
    defaultValue="parameter_file.as",
    allowedValues=["as"],
)
inputLigand = PluginVariable(
    name="Ligand SD file",
    id="input_ligand",
    description="The input ligand SD file.",
    type=VariableTypes.FILE,
    allowedValues=["sd","sdf"],
)

# ==========================#
# Variable outputs
# ==========================#
outputFile = PluginVariable(
    name="Output File",
    id="output_file",
    description="The output file with the ligand docked.",
    type=VariableTypes.FILE,
    defaultValue="output_dock",
)

##############################
#       Other variables      #
##############################
protocolPrmFile = PluginVariable(
    name="Protocol Parameter File",
    id="proto_prm_file",
    description="The docking protocol parameter file.",
    type=VariableTypes.FILE,
    defaultValue="dock.prm",
)
nRuns = PluginVariable(
    name="Number of runs",
    id="n_runs",
    description="Number of runs/ligand (default=10).",
    type=VariableTypes.INTEGER,
    defaultValue=10,
)
allH = PluginVariable(
    name="Hydrogens",
    id="all_h",
    description="Keep all hydrogens, read all hydrogens present (default=polar hydrogens only).",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
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
    if input_ligand is None:
        raise Exception("No ligand file provided.")
    if not os.path.exists(input_ligand):
        raise Exception("Ligand file does not exist.")
    
    cavity_file = block.inputs.get(cavityFile.id, None)
    if cavity_file is None:
        raise Exception("No cavity file provided.")
    if not os.path.exists(input_ligand):
        raise Exception("Cavity file does not exist.")

    out = "output_dock"
    path_output = os.path.dirname(input_ligand)
    path_cavity = os.path.dirname(cavity_file)

    if path_output != path_cavity:
        print("WARNING! The ligands to dock and the cavity file are not in the same directory.")
        print(f"The docking output file (.sd) will be saved in {path_output}.")

    output_file = os.path.join(path_output,block.outputs.get(outputFile.id, out))

    # rbcavity -was -d -r parameter_file.prm > parameter_file.log
    command = f"rbdock -i {input_ligand} -o {output_file} -r {input_PRMfile} "

    if block.variables.get("proto_prm_file", None) is not None:
        command += f"-p {block.variables.get('proto_prm_file')} "
    if block.variables.get("n_runs", 50) is not None:
        command += f"-n {block.variables.get('n_runs')} "
    if block.variables.get("all_h", True):
        command += "-allH "

    # Subprocess the command
    import subprocess

    completed_process = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Get the output and error
    output = completed_process.stdout
    error = completed_process.stderr

    # Save the output and error
    with open(f"{os.path.join(path_output,output_file)}.out", "w") as f:
        f.write(output)    
    with open(f"{os.path.join(path_output,output_file)}.err", "w") as f:
        f.write(error)

    # Set the output
    block.setOutput(outputFile.id, output_file)


rbDockBlock = PluginBlock(
    name="rDock",
    id="rbdock",
    description="Calculate the docking. (For local)",
    action=initialRbdock,
    variables=[
        protocolPrmFile,
        nRuns,
        allH,
    ],
    inputs=[inputPRMFile, cavityFile, inputLigand],
    outputs=[outputFile],
)
