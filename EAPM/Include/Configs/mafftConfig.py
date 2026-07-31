from HorusAPI import PluginConfig, PluginVariable, VariableTypes

mafftPathVariable = PluginVariable(
    id="mafft_path",
    name="MAFFT path",
    description="Path to the MAFFT executable",
    type=VariableTypes.FILE,
    defaultValue="MAFFT",
)


def checkMAFFTInstallation(block: PluginConfig):
    import os
    import shutil

    print("verifying MAFFT installation")

    # Get the path to the mafft executable
    mafftPath = block.variables.get(mafftPathVariable.id)

    # Check if the path is valid
    # Accept either an absolute path to the binary or a command available on the PATH
    if not (os.path.isfile(mafftPath) or shutil.which(mafftPath)):
        raise Exception(
            f"The MAFFT executable '{mafftPath}' was not found. "
            "Install MAFFT or set the correct path in the plugin configuration."
        )


# Create a plugin configuration for the mafft executable
mafftExecutableConfig = PluginConfig(
    name="MAFFT executable",
    description="Configure the path to the MAFFT executable for performing protein alignments",
    variables=[mafftPathVariable],
    action=checkMAFFTInstallation,
)
