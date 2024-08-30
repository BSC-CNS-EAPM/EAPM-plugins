"""
Module containing the rDock block for the EAPM plugin
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
    allowedValues=["sd", "sdf"],
)

# ==========================#
# Variable outputs
# ==========================#
outputFolder = PluginVariable(
    name="Docked Folder",
    id="output_file",
    description="The output folder with all the ligand docked in different files.",
    type=VariableTypes.FOLDER,
    defaultValue="dock_output",
)

##############################
#       Other variables      #
##############################
cpus = PluginVariable(
    name="CPUs",
    id="cpus",
    description="Number of cpus to be used in the docking.",
    type=VariableTypes.INTEGER,
    defaultValue=3,
)
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
    import subprocess

    def _splitInputLigands(input_ligand, cpus):
        """
        Split the input ligand file into multiple files.

        Parameters:
        ----------
        input_ligand: str
            The input ligand file
        cpus: int
            The number of CPUs to be used in the docking
        """

        num_cpus = os.cpu_count()
        if cpus > num_cpus:
            cpus = num_cpus
            print(
                f"WARNING! The number of CPUs requested is higher than the available CPUs. The number of CPUs will be set to the maximum available ({num_cpus})."
            )

        # setting paths
        path_output = os.path.dirname(input_ligand)
        path_split_ligands = os.path.join(path_output, "ligands")
        split_mols_file = os.path.join(path_output, ".splitMols.sh")
        split_mols_runner = os.path.join(path_output, ".split.sh")

        # Generating split file
        if not os.path.isfile(split_mols_file):
            with open(split_mols_file, "w") as filein:
                filein.writelines(
                    "#!/bin/bash\n"
                    "#Usage: splitMols.sh <input> #Nfiles <outputRoot>\n"
                    "fname=$1\n"
                    "nfiles=$2\n"
                    "output=$3\n"
                    "molnum=$(grep -c '$$$$' $fname)\n"
                    'echo " - $molnum molecules found"\n'
                    "echo \" - Dividing '$fname' into $nfiles files\"\n"
                    'echo " "\n'
                    "rawstep=`echo $molnum/$nfiles | bc -l`\n"
                    "let step=$molnum/$nfiles\n"
                    "if [ ! `echo $rawstep%1 | bc` == 0 ]; then\n"
                    "        let step=$step+1;\n"
                    "fi;\n"
                    "sdsplit -$step -o$output $1\n"
                )

        if not os.path.isdir(path_split_ligands):
            os.mkdir(path_split_ligands)

        # Generating splitted ligand files
        with open(split_mols_runner, "w") as fileout:
            fileout.writelines(
                f"bash .splitMols.sh {input_ligand} {cpus} {os.path.join(path_split_ligands,'split')}\n"
            )

        return split_mols_file, split_mols_runner

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

    # Splitting ligands
    output_folder_path = os.path.join(path_output, block.outputs.get(outputFolder.id, out))
    os.makedirs(output_folder_path, exist_ok=True)

    n_cpus = block.variables.get("cpus", 3)
    split_mols_file, split_mols_runner = _splitInputLigands(input_ligand, n_cpus)
    command_split = f"bash {os.path.join(path_output,'.split.sh')}"
    _ = subprocess.run(command_split, shell=True, capture_output=True, text=True)

    # Running rDock
    command_back = " "
    if block.variables.get("proto_prm_file", None) is not None:
        command_back += f"-p {block.variables.get('proto_prm_file')} "
    if block.variables.get("n_runs", 50) is not None:
        command_back += f"-n {block.variables.get('n_runs')} "
    if block.variables.get("all_h", True):
        command_back += "-allH "

    # Iterating over number of requested cpus
    processes = []
    for i in range(1, n_cpus + 1):
        command = f"rbdock -i {os.path.join(path_output, f'ligands/split{i}.sd')} -o {os.path.join(output_folder_path, f'split{i}_out')} -r {input_PRMfile}"
        full_command = command + command_back

        # Start the process and store the Popen object
        process = subprocess.Popen(
            full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        processes.append((i, process))

    with open(os.path.join(output_folder_path, ".docking.info"), "w") as f:
        f.write(full_command)

    # Wait for all processes to complete and collect their outputs
    for i, process in processes:
        output, error = process.communicate()

        # Save the output and error
        with open(os.path.join(output_folder_path, f"dock_{i}.out"), "w") as f:
            f.write(output.decode())
        with open(os.path.join(output_folder_path, f"dock_{i}.err"), "w") as f:
            f.write(error.decode())

    # Removing the split files
    os.remove(split_mols_file)
    os.remove(split_mols_runner)

    # Set the output
    block.setOutput(outputFolder.id, output_folder_path)


rbDockBlock = PluginBlock(
    name="rDock",
    id="rbdock",
    description="Perform a docking with rDock. (For local)",
    action=initialRbdock,
    variables=[
        cpus,
        protocolPrmFile,
        nRuns,
        allH,
    ],
    inputs=[inputPRMFile, cavityFile, inputLigand],
    outputs=[outputFolder],
)
