import datetime
import os
import subprocess

from HorusAPI import PluginBlock, PluginVariable, VariableList, VariableTypes

# ==========================#
# Variable inputs
# ==========================#
dbPath = PluginVariable(
    name="Database path",
    id="db_path",
    description="The database where the search will be done.",
    type=VariableTypes.FILE,
    defaultValue="nr.fa",
    allowedValues=["fa"],
)
containerName = PluginVariable(
    name="Container name",
    id="container_name",
    description="The container name to use.",
    type=VariableTypes.STRING,
    defaultValue="bsceapm/ahatool:2.2",
)
inputFasta = PluginVariable(
    name="Input fasta",
    id="input_fasta",
    description="The input fasta file.",
    type=VariableTypes.FILE,
    defaultValue=None,
    allowedValues=["fasta"],
)

# ==========================#
# Variable outputs
# ==========================#
outputAhatool = PluginVariable(
    name="Ahatool output",
    id="path",
    description="The folder containing the results.",
    type=VariableTypes.FOLDER,
)

##############################
#       Other variables      #
##############################
prefixVar = PluginVariable(
    name="Prefix",
    id="prefix",
    description="The prefix the tool will use for produced files.",
    type=VariableTypes.STRING,
    defaultValue=None,
)
startVar = PluginVariable(
    name="Start",
    id="start",
    description="Start of execution (search or build).",
    type=VariableTypes.STRING,
    defaultValue=None,
)
evalVar = PluginVariable(
    name="E-value",
    id="e_value",
    description="E-value (recommended: 1e-10).",
    type=VariableTypes.FLOAT,
    defaultValue=None,
)
threadsVar = PluginVariable(
    name="Threads",
    id="threads",
    description="Number of threads. Options (1, 2, 4).",
    type=VariableTypes.INTEGER,
    defaultValue=None,
)


def initialAction(block: PluginBlock):

    container_name = block.inputs.get("container_name", "bsceapm/ahatool:2.2")
    input_fasta = block.inputs.get("input_fasta", None)
    fasta = os.path.basename(input_fasta)
    fasta_root, _ = os.path.split(input_fasta)
    database_path = block.inputs.get("db_path", None)
    database_name = database_path.split("/")[-1]
    database_root, _ = os.path.split(database_path)

    # Test input fasta
    if input_fasta is None:
        raise Exception("No input fasta provided")
    if not os.path.exists(input_fasta):
        raise Exception(f"The input fasta file does not exist: {input_fasta}")

    # Test docker installation
    try:
        subprocess.check_output(["docker", "--version"])
    except subprocess.CalledProcessError:
        raise Exception("Docker is not installed")

    # Test ahatool docker
    try:
        output = subprocess.check_output(["docker", "images", "-q", container_name])
    except subprocess.CalledProcessError:
        raise Exception("Failed to list Docker images")

    if not output:
        raise Exception(
            f"Docker image not found: {container_name}. Please pull the image first. \n\t !!!!For academic use only!!!!"
        )

    start = block.variables.get("start", None)
    prefix = block.variables.get("prefix", None)
    eval = block.variables.get("e_value", None)
    threads = block.variables.get("threads", None)

    # docker run -v /home/albertcs/GitHub/NBD/test/ahatool:/home/projects bsceapm/ahatool:2.2 /bin/bash ../AHATools/AHATool.sh -i /home/albertcs/GitHub/NBD/test/ahatool/test.fasta -d nr_test.fa -p ahatool_output_1715175010.575412

    if prefix is None:
        output_name = "ahatool_output_" + str(datetime.datetime.now().timestamp())
    else:
        output_name = prefix + "_" + str(datetime.datetime.now().timestamp())

    command = f"docker run -v {fasta_root}:/home/projects -v {database_root}:/home/database {container_name} /bin/bash ../AHATool/AHATool.sh -i {fasta} -d {database_name} "
    if start is not None:
        command += f"-s {start} "
    if output_name is not None:
        command += f"-p {output_name} "
    if eval is not None:
        command += f"-e {eval} "
    if threads is not None:
        command += f"-t {threads} "

    # Run with subprocess
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(f"Failed to run command: {result.stderr}")
        output = result.stdout
    except subprocess.CalledProcessError as e:
        raise Exception(f"Failed to run command: {e}")
    block.setOutput(outputAhatool.id, output_name)


ahatoolBlock = PluginBlock(
    name="AHATool",
    action=initialAction,
    description="Iteratively search a protein sequence against a protein database",
    inputs=[inputFasta, dbPath, containerName],
    variables=[prefixVar, startVar, evalVar, threadsVar],
    outputs=[outputAhatool],
)
