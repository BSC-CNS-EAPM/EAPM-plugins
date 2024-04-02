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

    print("verifying HMMER installation")

    # Get the path to the mafft executable
    hmmerPath = block.variables.get("HMMER_path")

    # Check if the path is valid
    if not os.path.isfile(hmmerPath):
        raise Exception("The HMMER executable path is not valid")


# Create a plugin configuration for the mafft executable
hmmerExecutableConfig = PluginConfig(
    name="HMMER executable",
    description="Configure the path to the HMMER executable for performing protein alignments",
    variables=[hmmerPathVariable],
    action=checkHmmerInstallation,
)
