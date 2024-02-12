import os
import shutil
import datetime
import subprocess
import typing

from HorusAPI import SlurmBlock, PluginBlock, VariableList, PluginVariable, VariableTypes

localIPs = {"cactus": "84.88.51.217", "blossom": "84.88.51.250", "bubbles": "84.88.51.219"}


def setup_bsc_calculations_based_on_horus_remote(
    remote_name, remote_host: str, jobs, partition, scriptName, cpus, job_name, program, modulePurge
):
    import bsc_calculations

    cluster = "local"

    if remote_name != "local":
        cluster = remote_host

    if remote_host in localIPs.values():
        cluster = "powerpuff"

    # If we are working with pele, only marenostrum and nord3 are allowed
    if program == "pele":
        if cluster not in ["mn1.bsc.es", "mn2.bsc.es", "mn3.bsc.es", "nord3.bsc.es"]:
            raise Exception("Pele can only be run on Marenostrum or Nord3")

        if cluster == "nord3.bsc.es":
            bsc_calculations.nord3.setUpPELEForNord3(
                jobs,
                partition=partition,
                cpus=cpus,
                general_script=scriptName,
                scripts_folder=scriptName + "_scripts",
            )
        elif "mn" in cluster:
            bsc_calculations.marenostrum.setUpPELEForMarenostrum(
                jobs,
                partition=partition,
                cpus=cpus,
                general_script=scriptName,
                scripts_folder=scriptName + "_scripts",
            )

        return cluster

    ## Define cluster
    # cte_power
    if cluster == "plogin1.bsc.es":
        bsc_calculations.cte_power.jobArrays(
            jobs,
            job_name=job_name,
            partition=partition,
            program=program,
            script_name=scriptName,
            gpus=cpus,
            module_purge=modulePurge,
        )
    # marenostrum
    elif "mn" in cluster:
        bsc_calculations.marenostrum.jobArrays(
            jobs,
            job_name=job_name,
            partition=partition,
            program=program,
            script_name=scriptName,
            cpus=cpus,
            module_purge=modulePurge,
        )
    # minotauro
    elif cluster == "mt1.bsc.es":
        bsc_calculations.minotauro.jobArrays(
            jobs,
            job_name=job_name,
            partition=partition,
            program=program,
            script_name=scriptName,
            gpus=cpus,
            module_purge=modulePurge,
        )
    # powerpuff
    elif cluster == "powerpuff":
        print("Generating powerpuff girls jobs...")
        bsc_calculations.local.parallel(
            jobs,
            cpus=min(40, len(jobs)),
            script_name=scriptName,
        )
    # local
    elif cluster == "local":
        print("Generating local jobs...")
        bsc_calculations.local.parallel(
            jobs,
            cpus=min(40, len(jobs)),
            script_name=scriptName,
        )
    else:
        raise Exception("Cluster not supported.")

    return cluster


HOOK_SCRIPT = """
for script in calculation_script.sh_?; do
    sh "$script" > "${script%.*}.out" 2> "${script%.*}.err" &
    exit_code=$?
done

# Wait for all background processes to finish
wait

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


def launchCalculationAction(
    block: SlurmBlock,
    jobs: typing.List[str],
    program: str,
    uploadFolders: typing.Optional[typing.List[str]] = None,
    modulePurge: typing.Optional[bool] = False,
):
    if jobs is None:
        raise Exception("No jobs selected")

    partition = block.variables.get("partition")
    cpus = block.variables.get("cpus")
    simulationName = block.variables.get("folder_name")
    scriptName = block.variables.get("script_name", "calculation_script.sh")

    if simulationName is None:
        simulationName = block.flow.name.lower().replace(" ", "_")

    block.extraData["simulationName"] = simulationName

    print(f"Launching BSC calculation with {cpus} CPUs")

    cluster = setup_bsc_calculations_based_on_horus_remote(
        block.remote.name.lower(),
        block.remote.host,
        jobs,
        partition,
        scriptName,
        cpus,
        simulationName,
        program,
        modulePurge
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
        savedID_and_date = block.flow.savedID + "_" + str(datetime.datetime.now().timestamp())
        simRemoteDir = os.path.join(block.remote.workDir, savedID_and_date)
        block.extraData["remoteDir"] = simRemoteDir
        block.remote.remoteCommand(f"mkdir -p -v {simRemoteDir}")

        print(f"Created simulation folder in the remote at {simRemoteDir}")
        print("Sending data to the remote...")

        # Check if in the input, scpefic folders to upload are specified
        # If so, upload them
        if uploadFolders is not None:
            for file in uploadFolders:
                finalPath = block.remote.sendData(file, simRemoteDir)
            block.extraData["uploadedFolder"] = False
        else:
            # Send the whole folder to the remote
            simRemoteDir = block.remote.sendData(os.getcwd(), simRemoteDir)
            block.extraData["uploadedFolder"] = True

        block.extraData["remoteContainer"] = simRemoteDir
        # base_folder = os.path.basename(os.getcwd())

        # # Move the contents of the sent folder to its parent
        # # This is done because the folder is sent as a subfolder
        # command = f"command: mv {simRemoteDir}/{base_folder} {simRemoteDir}"
        # block.remote.remoteCommand(f"mv {simRemoteDir}/{base_folder} {simRemoteDir}")

        # # Remove the sent folder
        # block.remote.remoteCommand(f"rm -rf {simRemoteDir}/{base_folder}")

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
            print(f"Submitting the job to the remote... {scriptPath}")
            if program == "pele":
                with block.remote.cd(simRemoteDir):
                    for jobScript in os.listdir(scriptName + "_scripts"):
                        if jobScript.endswith(".sh"):
                            jobID = block.remote.submitJob(
                                scriptPath + "_scripts/" + jobScript, changeDir=False
                            )
                            print("Submitted job with ID: ", jobID)
                print("Waiting for the jobs to finish...")
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

        # Create the folder
        os.makedirs(folderDestinationOverride)

        final_path = block.remote.getData(simRemoteDir, folderDestinationOverride)

        # If we sent the whole folder, the results are in a subfolder
        # Move them to the parent folder
        if block.extraData.get("uploadedFolder", False):
            print("Uploaded folder, moving results to parent folder")
            final_path = os.path.join(final_path, os.path.basename(currentFolder))

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

        remove_remote_folder_on_finish = block.variables.get("remove_folder_on_finish", True)
        # Remove the remote folder
        if remove_remote_folder_on_finish:
            print(f"Removing remote folder {remoteContainer}")
            block.remote.remoteCommand(f"rm -rf {remoteContainer}")
    else:
        final_path = os.path.join(os.getcwd())
        print("Calculation finished, results are in the folder: ", final_path)

    return final_path


# Other variables
simulationNameVariable = PluginVariable(
    name="Simulation name",
    id="folder_name",
    description="Name of the simulation folder. By default it will be the same as the flow name.",
    type=VariableTypes.STRING,
    category="Slurm configuration",
)
scriptNameVariable = PluginVariable(
    name="Script name",
    id="script_name",
    description="Name of the script.",
    type=VariableTypes.STRING,
    defaultValue="calculation_script.sh",
    category="Slurm configuration",
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

removeFolderOnFinishVariable = PluginVariable(
    name="Remove remote folder on finish",
    id="remove_folder_on_finish",
    description="Deletes the calculation folder on the remote on finish.",
    type=VariableTypes.BOOLEAN,
    defaultValue=True,
    category="Remote",
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

BSC_JOB_VARIABLES = [
    simulationNameVariable,
    scriptNameVariable,
    partitionVariable,
    cpusVariable,
    environmentList,
    removeFolderOnFinishVariable,
]
