import os
import shutil
import subprocess

from HorusAPI import SlurmBlock, VariableTypes, PluginVariable, VariableList

# Input variables
inputJobsVariable = PluginVariable(
    id="input_jobs",
    name="Jobs",
    description="Jobs for the BSC calculations library",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_jobs"],
)

# Other variables
simulationNameVariable = PluginVariable(
    name="Simulation name",
    id="folder_name",
    description="Name of the simulation folder. By default it will be the same as the flow name.",
    type=VariableTypes.STRING,
)
scriptNameVariable = PluginVariable(
    name="Script name",
    id="script_name",
    description="Name of the script.",
    type=VariableTypes.STRING,
    defaultValue="calculation_script.sh",
)

partitionVariable = PluginVariable(
    name="Partition",
    id="partition",
    description="Partition where to lunch.",
    type=VariableTypes.STRING_LIST,
    defaultValue="bsc_ls",
    allowedValues=["bsc_ls", "debug"],
    category="Slurm configuration",
)

cpusVariable = PluginVariable(
    name="CPUs",
    id="cpus",
    description="Number of CPUs to use.",
    type=VariableTypes.INTEGER,
    defaultValue=1,
    category="Slurm configuration",
)

# Advanced variables
environmentKeyVariable = PluginVariable(
    name="Environment",
    id="environment_key",
    description="Environment key",
    type=VariableTypes.STRING,
    category="Environment",
)

environmentValueVariable = PluginVariable(
    name="Value",
    id="environment_value",
    description="Environment value",
    type=VariableTypes.STRING,
    category="Environment",
)

environmentList = VariableList(
    id="environment_list",
    name="Environment variables",
    description="Environment variables to set during the remote connection.",
    prototypes=[environmentKeyVariable, environmentValueVariable],
    category="Environment",
)

# Output variables
outputSimulationResultsVariable = PluginVariable(
    id="output_simulation_results",
    name="Results",
    description="Output simulation results",
    type=VariableTypes.CUSTOM,
    allowedValues=["bsc_results"],
)

outputResultsFolderVariable = PluginVariable(
    id="output_results_folder",
    name="Results folder",
    description="Output results folder",
    type=VariableTypes.FOLDER,
)


def launchCalculationAction(block: SlurmBlock):
    jobs_input = block.inputs.get("input_jobs")

    program = jobs_input["program"]
    jobs = jobs_input["jobs"]

    if jobs is None:
        raise Exception("No jobs selected")

    partition = block.variables.get("partition")
    cpus = block.variables.get("cpus")
    simulationName = block.variables.get("folder_name")
    scriptName = block.variables.get("script_name", "calculation_script.sh")

    if simulationName is None:
        simulationName = block.flow.name.lower().replace(" ", "_")

    block.extraData["simulationName"] = simulationName

    from Utils.bsc_calculations_cluster import setup_bsc_calculations_based_on_horus_remote

    print("Launching BSC calculations")

    cluster = setup_bsc_calculations_based_on_horus_remote(
        block.remote.name.lower(),
        block.remote.host,
        jobs,
        partition,
        scriptName,
        cpus,
        simulationName,
        program,
    )

    # Read the environment variables
    environmentValues = block.variables.get("environment_list", [])
    environmentListValues = {}
    if environmentValues is not None:
        for env in environmentValues:
            environmentListValues[env["environment_key"]] = env["environment_value"]

    # Rewrite the main script to add the environment variables
    # and allow for waiting for the jobs to finish
    # This is only necessary for powerpuff and local
    if cluster == "powerpuff" or cluster == "local":
        with open(scriptName, "w") as f:
            f.write("#!/bin/sh\n")

            for key, value in environmentListValues.items():
                f.write(f"export {key}={value}\n")

            f.write(HOOK_SCRIPT)

    if cluster != "local":
        simRemoteDir = os.path.join(block.remote.workDir, block.flow.savedID)
        block.extraData["remoteDir"] = simRemoteDir
        block.remote.remoteCommand(f"mkdir -p -v {simRemoteDir}")
        block.extraData["remoteContainer"] = simRemoteDir

        print(f"Created simulation folder in the remote at {simRemoteDir}")
        print("Sending data to the remote...")

        # Check if in the input, scpefic folders to upload are specified
        # If so, upload them
        toUpload = jobs_input.get("send_data", None)

        if toUpload is not None:
            for file in toUpload:
                block.remote.sendData(file, simRemoteDir)
        else:
            # Send the whole folder to the remote
            block.remote.sendData(os.getcwd(), simRemoteDir)

            base_folder = os.path.basename(os.getcwd())

            # Move the contents of the sent folder to its parent
            # This is done because the folder is sent as a subfolder
            block.remote.remoteCommand(f"mv {simRemoteDir}/{base_folder} {simRemoteDir}")

            # Remove the sent folder
            block.remote.remoteCommand(f"rm -rf {simRemoteDir}/{base_folder}")

        # Upload the commands
        for file in os.listdir("."):
            if file.startswith(scriptName):
                block.remote.sendData(file, simRemoteDir)

        # Upload the script
        scriptPath = block.remote.sendData(scriptName, simRemoteDir)

        print("Data sent to the remote.")

        print("Running the simulation...")

        # Run the simulation
        if cluster == "powerpuff":
            # The powerpuff cluster doesn't have Slurm, so we need to run the script manually & load the Schrodinger module
            schrodingerPath = block.remote.remoteCommand("echo $SCHRODINGER")
            command = f"export={schrodingerPath} cd {simRemoteDir} && bash {scriptName}"

            block.remote.remoteCommand(command)

        else:
            jobID = block.remote.submitJob(scriptPath)

            print(f"Simulation running with job ID {jobID}. Waiting for it to finish...")
    # * Local
    else:
        print("Running the simulation locally...")

        oldEnv = os.environ.copy()

        for key, value in environmentListValues.items():
            os.environ[key] = value

        # Run the simulation
        try:
            with subprocess.Popen(
                ["sh", scriptName],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            ) as p:
                print(f"Simulation running with PID {p.pid}. Waiting for it to finish...")

                if p.stdout is None:
                    raise Exception("No stdout produced by the process")

                for line in p.stdout:
                    strippedOut = line.decode("utf-8").strip()
                    if strippedOut != "":
                        print(strippedOut)

                # Print the error
                strippedErr = ""
                if p.stderr:
                    for line in p.stderr:
                        strippedErr = line.decode("utf-8").strip()
                        if strippedErr != "":
                            print(strippedErr)

                # Wait for the process to finish
                p.wait()

                if p.returncode != 0:
                    raise Exception(strippedErr)
        finally:
            os.environ = oldEnv


def downloadResultsAction(block: SlurmBlock):
    """
    Final action of the block. It downloads the results from the remote.

    Args:
        block (SlurmBlock): The block to run the action on.
    """

    if block.remote.name != "Local":
        cluster = block.remote.host
    else:
        cluster = "local"

    if cluster != "local":
        simRemoteDir = block.extraData["remoteDir"]

        print("Calculation finished, downloading results...")

        currentFolder = os.getcwd()
        folderDestinationOverride = os.path.join(currentFolder, "tmp_download")

        if os.path.exists(folderDestinationOverride):
            shutil.rmtree(folderDestinationOverride)

        final_path = block.remote.getData(simRemoteDir, folderDestinationOverride)

        # Move the contents of the downloaded folder to its parent
        # This is done because the folder is downloaded as a subfolder
        for file in os.listdir(final_path):
            current_path = os.path.join(final_path, file)
            new_path = os.path.join(currentFolder, file)

            if os.path.exists(new_path):
                if os.path.isdir(new_path):
                    shutil.rmtree(new_path)
                else:
                    os.remove(new_path)

            shutil.move(current_path, new_path)

        # Remove the downloaded folder
        shutil.rmtree(folderDestinationOverride)

        final_path = currentFolder

        print(f"Results downloaded to {final_path}")

        remoteContainer = block.extraData["remoteContainer"]

        # Remove the remote folder
        print(f"Removing remote folder {remoteContainer}")
        block.remote.remoteCommand(f"rm -rf {remoteContainer}")
    else:
        final_path = os.path.join(os.getcwd())
        print("Calculation finished, results are in the folder: ", final_path)

    block.setOutput("output_results_folder", final_path)

    # Set as the output results the input configuration
    # This is very specific to EAPM blocks
    jobs_input = block.inputs.get("input_jobs")

    jobs_input = jobs_input.get("results_data")

    block.setOutput("output_simulation_results", jobs_input)


bscCalculationsBlock = SlurmBlock(
    name="BSC calculations",
    description="Launch BSC calculations jobs",
    inputs=[
        inputJobsVariable,
    ],
    variables=[
        partitionVariable,
        cpusVariable,
        simulationNameVariable,
        scriptNameVariable,
        environmentList,
    ],
    outputs=[outputSimulationResultsVariable, outputResultsFolderVariable],
    initialAction=launchCalculationAction,
    finalAction=downloadResultsAction,
)


HOOK_SCRIPT = """
for script in calculation_script.sh_?; do
    sh "$script" > "${script%.*}.out" 2> "${script%.*}.err"
    exit_code=$?
    
    # Wait for the last background process to finish
    wait $!
done

if [ $exit_code -ne 0 ]; then
    echo "Error: Script $script failed with exit code $exit_code" >&2
    exit 1
fi

# Check if the .err file is empty in order to determine
# if the script ran successfully
if [ -s "${script%.*}.err" ]; then
    echo "Error: Script $script failed with errors:" >&2
    cat "${script%.*}.err" >&2
    exit 1
fi


echo "All scripts completed successfully."

"""
