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

    print("verifying MAFFT installation")

    # Get the path to the mafft executable
    mafftPath = block.variables.get("MAFFT_path")

    # Check if the path is valid
    if not os.path.isfile(mafftPath):
        raise Exception("The MAFFT executable path is not valid")


# Create a plugin configuration for the mafft executable
mafftExecutableConfig = PluginConfig(
    name="MAFFT executable",
    description="Configure the path to the MAFFT executable for performing protein alignments",
    variables=[mafftPathVariable],
    action=checkMAFFTInstallation,
)
