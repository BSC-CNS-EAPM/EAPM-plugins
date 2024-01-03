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

    from Blocks.BSCCalculations import bscCalculationsBlock

    eapmPlugin.addBlock(bscCalculationsBlock)

    from Blocks.AlphaFoldEAPM import alphafoldBlock  # type: ignore

    # Add the block to the plugin
    eapmPlugin.addBlock(alphafoldBlock)

    from Blocks.PrepWizardEAPM import prepWizardBlock  # type: ignore

    eapmPlugin.addBlock(prepWizardBlock)

    from Blocks.AlignPdbEAPM import alignBlock  # type: ignore

    eapmPlugin.addBlock(alignBlock)

    from Blocks.SetupDockingGrid import setupDockingGrid  # type: ignore

    eapmPlugin.addBlock(setupDockingGrid)

    from Blocks.AnalyseGlideDocking import analyseGlideDocking

    eapmPlugin.addBlock(analyseGlideDocking)

    from Blocks.PDBToMAE import convertPDBToMAEBlock

    eapmPlugin.addBlock(convertPDBToMAEBlock)

    from Blocks.TrimAlphafoldModels import trimAlphaFoldModelsBlock

    eapmPlugin.addBlock(trimAlphaFoldModelsBlock)

    from Blocks.SetupGlide import setupGlideBlock

    eapmPlugin.addBlock(setupGlideBlock)

    # Return the plugin
    return eapmPlugin


plugin = createPlugin()
