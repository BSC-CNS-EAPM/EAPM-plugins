from HorusAPI import PluginConfig, PluginVariable, VariableTypes

pyrosettaPythonVariable = PluginVariable(
    id="pyrosetta_python",
    name="PyRosetta python",
    description="Path to the python interpreter that has PyRosetta installed. "
    "It is used to run the Rosetta analysis scripts.",
    type=VariableTypes.FILE,
    defaultValue="python",
)


def checkPyrosettaInstallation(block: PluginConfig):
    import os
    import shutil
    import subprocess

    print("verifying PyRosetta installation")

    pythonPath = block.variables.get(pyrosettaPythonVariable.id)

    # Accept either an absolute path to the interpreter or a command on the PATH
    resolved = pythonPath if os.path.isfile(pythonPath) else shutil.which(pythonPath)

    if not resolved:
        raise Exception(
            f"The python interpreter '{pythonPath}' was not found. "
            "Set the path to a python with PyRosetta in the plugin configuration."
        )

    # Verify that the interpreter can actually import pyrosetta
    result = subprocess.run(
        [resolved, "-c", "import pyrosetta"],
        capture_output=True,
        check=False,
    )

    if result.returncode != 0:
        raise Exception(
            f"The python interpreter '{resolved}' cannot import PyRosetta. "
            "Install PyRosetta in that environment or point the configuration "
            "to the environment where it is installed."
        )


# Create a plugin configuration for the python with PyRosetta
pyrosettaExecutableConfig = PluginConfig(
    name="PyRosetta python",
    description="Configure the python interpreter used to run the Rosetta analysis scripts",
    variables=[pyrosettaPythonVariable],
    action=checkPyrosettaInstallation,
)
