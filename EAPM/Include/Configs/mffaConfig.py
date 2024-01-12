from HorusAPI import PluginConfig, PluginVariable, VariableTypes

mffaPathVariable = PluginVariable(
    id="mffa_path",
    name="MFFA path",
    description="Path to the MFFA executable",
    type=VariableTypes.FILE,
    defaultValue="mffa",
)


def checkMFFAInstallation(block: PluginConfig):
    import os

    print("verifying mffa installation")

    # Get the path to the mffa executable
    mffaPath = block.variables.get("mffa_path")

    # Check if the path is valid
    if not os.path.isfile(mffaPath):
        raise Exception("The MFFA executable path is not valid")


# Create a plugin configuration for the mffa executable
mffaExecutableConfig = PluginConfig(
    name="MFFA executable",
    description="Configure the path to the MFFA executable for performing protein alignments",
    variables=[mffaPathVariable],
    action=checkMFFAInstallation,
)
