"""
Plugin configuration for the MMseqs2 executable used by the clustering block.
"""

from HorusAPI import PluginConfig, PluginVariable, VariableTypes

mmseqsPathVariable = PluginVariable(
    id="mmseqs_path",
    name="MMseqs2 path",
    description="Path to the MMseqs2 executable (or the command name if it is on the PATH)",
    type=VariableTypes.FILE,
    defaultValue="mmseqs",
)


def checkMMseqsInstallation(block: PluginConfig):
    import os
    import shutil

    print("Verifying MMseqs2 installation")

    mmseqs_path = block.variables.get(mmseqsPathVariable.id, "mmseqs")

    # Accept either an absolute path to the binary or a command available on the PATH
    if not (os.path.isfile(mmseqs_path) or shutil.which(mmseqs_path)):
        raise Exception(
            f"The MMseqs2 executable '{mmseqs_path}' was not found. "
            "Install MMseqs2 or set the correct path in the plugin configuration."
        )


# Create a plugin configuration for the mmseqs executable
mmseqsExecutableConfig = PluginConfig(
    name="MMseqs2 executable",
    description="Configure the path to the MMseqs2 executable for sequence clustering",
    variables=[mmseqsPathVariable],
    action=checkMMseqsInstallation,
)
