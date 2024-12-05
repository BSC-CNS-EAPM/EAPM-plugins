"""
Module containing the rDock block for the rDock plugin
"""

from HorusAPI import SlurmBlock, PluginVariable, VariableTypes, PluginConfig
import os
import subprocess
import shutil
import datetime

rDockPath_variable = PluginVariable(
    name="rDock Path",
    id="rdock_path",
    description="The path to the rDock installation.",
    type=VariableTypes.FILE,
)

rDock_path = PluginConfig(
    name="rDock Path",
    id="rdock_path",
    description="The path to the rDock installation.",
    variables=[rDockPath_variable],
)

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
    description="The input ligand or ligands in SD file.",
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
    description="Number of docking poses per ligand (default=10).",
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

SLURM_DEFAULT_PREAMBLE = """#!/bin/bash
#SBATCH -J rdock
#SBATCH --output=rdock_%j.out
#SBATCH --error=rdock_%j.err
#SBATCH --ntasks=1
#SBATCH --time=30:00
#SBATCH --cpus-per-task 1
#SBATCH --array=1-%n_runs%

# Load necessary modules as needed

# The rDock block will place here the execution command
%execution_command%
"""

slurmScript = PluginVariable(
    name="Slurm script",
    id="slurm_script",
    description="The slurm script to be executed.",
    type=VariableTypes.CODE,
    defaultValue=SLURM_DEFAULT_PREAMBLE,
    allowedValues=["shell"],
)


def _localRunRDock(
    rdock_executable_path, input_PRMfile, n_cpus, path_output, output_folder_path, command_back
):

    # Iterating over number of requested cpus
    processes = []
    for i in range(1, n_cpus + 1):
        command = f"{rdock_executable_path} -i {os.path.join(path_output, f'ligands/split{i}.sd')} -o {os.path.join(output_folder_path, f'split{i}_out')} -r {input_PRMfile}"
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


# Align action block
def initialRbdock(block: SlurmBlock):

    import os
    import subprocess
    import stat

    def _splitInputLigands(input_ligand, cpus, path_split_ligands):
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
        path_output = os.getcwd()
        split_mols_file = os.path.join(path_output, ".splitMols.sh")

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

        command = f"bash {split_mols_file} {input_ligand} {cpus} {os.path.join(path_split_ligands,'split')}\n"
        os.chmod(split_mols_file, os.stat(split_mols_file).st_mode | stat.S_IEXEC)

        return command

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
    path_output = os.getcwd()

    # Splitting ligands
    output_folder_path = os.path.join(path_output, out)
    os.makedirs(output_folder_path, exist_ok=True)

    n_cpus = block.variables.get("cpus", 3)
    path_split_ligands = os.path.join(path_output, "ligands")
    command = _splitInputLigands(input_ligand, n_cpus, path_split_ligands)
    subprocess.run(
        command,  # Pass the command as a string
        shell=True,  # Enable shell mode
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Running rDock
    command_back = " "
    if block.variables.get("proto_prm_file", None) is not None:
        command_back += f"-p {block.variables.get('proto_prm_file')} "
    if block.variables.get("n_runs", 50) is not None:
        command_back += f"-n {block.variables.get('n_runs')} "
    if block.variables.get("all_h", True):
        command_back += "-allH "

    rdock_executable_path = block.config.get("rdock_path", None)

    if block.remote.isLocal:
        _localRunRDock(
            rdock_executable_path,
            input_PRMfile,
            n_cpus,
            path_output,
            output_folder_path,
            command_back,
        )
    else:
        sending_folder = f"rdock_sandbox_{datetime.datetime.now().timestamp()}"
        if os.path.exists(sending_folder):
            shutil.rmtree(sending_folder)
        os.makedirs(sending_folder, exist_ok=True)

        # Copying input files to the sandbox folder
        shutil.copy(input_PRMfile, sending_folder)
        shutil.copytree(path_split_ligands, sending_folder + "/ligands")
        shutil.copy(cavity_file, sending_folder)

        print("Sending data to remote server...")
        print(sending_folder + ".tar.gz", block.remote.workDir)
        remote_folder = block.remote.sendData(sending_folder, block.remote.workDir)
        print("Data sent.")
        input_remote_PRMfile = os.path.join(remote_folder, os.path.basename(input_PRMfile))
        output_remote_folder = os.path.join(remote_folder, os.path.basename(output_folder_path))

        block.remote.remoteCommand(f"mkdir -p {output_remote_folder}")

        # Generate slurm script
        slurm_value = block.variables.get(slurmScript.id, SLURM_DEFAULT_PREAMBLE)

        slurm_command = ""
        for i in range(1, n_cpus + 1):
            slurm_command += f"if [[ $SLURM_ARRAY_TASK_ID = {i} ]]; then\n"

            command = f"{rdock_executable_path} -i {os.path.join(remote_folder, f'ligands/split{i}.sd')} -o {os.path.join(output_remote_folder, f'split{i}_out')} -r {input_remote_PRMfile}"
            slurm_command += command + command_back + "\n"

            slurm_command += "fi\n"

        slurm_value = slurm_value.replace("%execution_command%", slurm_command)

        # Write the slurm script
        with open("rdock.slurm", "w") as f:
            f.write(slurm_value)

        slurm_remote = block.remote.sendData("rdock.slurm", remote_folder)
        jobId = block.remote.submitJob(slurm_remote)

        print(f"Job submitted with ID: {jobId}")

        block.extraData["remote_folder"] = remote_folder


def download_results(block: SlurmBlock):

    output_folder_path = "output_dock"

    if not block.remote.isLocal:
        remote_folder = block.extraData.get("remote_folder", None)
        if remote_folder is None:
            raise Exception("No remote folder found.")

        block.remote.getData(remote_folder, output_folder_path)

    # Set the output
    block.setOutput(outputFolder.id, output_folder_path)


rbDockBlock = SlurmBlock(
    name="rDock",
    id="rbdock",
    description="Perform a docking with rDock. (For local)",
    initialAction=initialRbdock,
    finalAction=download_results,
    variables=[
        cpus,
        protocolPrmFile,
        nRuns,
        allH,
        slurmScript,
    ],
    inputs=[inputPRMFile, cavityFile, inputLigand],
    outputs=[outputFolder],
)
