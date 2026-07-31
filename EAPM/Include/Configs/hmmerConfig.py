from HorusAPI import PluginConfig, PluginVariable, VariableTypes

hmmerPathVariable = PluginVariable(
    id="hmmer_path",
    name="HMMER path",
    description="Path to the HMMER executable",
    type=VariableTypes.FILE,
    defaultValue="HMMER",
)


def checkHmmerInstallation(block: PluginConfig):
    import os
    import shutil

    print("verifying HMMER installation")

    # Get the path to the mafft executable
    hmmerPath = block.variables.get(hmmerPathVariable.id)

    # Check if the path is valid
    # Accept either an absolute path to the binary or a command available on the PATH
    if not (os.path.isfile(hmmerPath) or shutil.which(hmmerPath)):
        raise Exception(
            f"The HMMER executable '{hmmerPath}' was not found. "
            "Install HMMER or set the correct path in the plugin configuration."
        )


# Create a plugin configuration for the mafft executable
hmmerExecutableConfig = PluginConfig(
    name="HMMER executable",
    description="Configure the path to the HMMER executable for performing protein alignments",
    variables=[hmmerPathVariable],
    action=checkHmmerInstallation,
)
