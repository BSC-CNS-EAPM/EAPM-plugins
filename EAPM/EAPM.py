"""
Entry point for the EAPM plugin
"""

from HorusAPI import Plugin


def createPlugin():
    """
    Generates the EAPM plugin and returns the instance
    """
    # ========== Plugin Definition ========== #

    eapmPlugin = Plugin(id="EAPM")

    # ========== Blocks ========== #

    from Blocks.AlphaFoldEAPM import alphafoldBlock  # type: ignore

    # Add the block to the plugin
    eapmPlugin.addBlock(alphafoldBlock)

    from Blocks.PrepWizardEAPM import prepWizardBlock  # type: ignore

    eapmPlugin.addBlock(prepWizardBlock)

    # Return the plugin
    return eapmPlugin


plugin = createPlugin()
