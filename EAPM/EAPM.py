"""
Entry point for the EAPM plugin
"""

from HorusAPI import Plugin, PluginConfig


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

    from Blocks.AlignPdbEAPM import alignBlock  # type: ignore

    eapmPlugin.addBlock(alignBlock)

    from Blocks.SetupDockingGrid import setupDockingGrid  # type: ignore

    eapmPlugin.addBlock(setupDockingGrid)

    from Blocks.AnalyseGlideDocking import analyseGlideDocking  # type: ignore

    eapmPlugin.addBlock(analyseGlideDocking)

    from Blocks.PDBToMAE import convertPDBToMAEBlock  # type: ignore

    eapmPlugin.addBlock(convertPDBToMAEBlock)

    from Blocks.TrimAlphafoldModels import trimAlphaFoldModelsBlock  # type: ignore

    eapmPlugin.addBlock(trimAlphaFoldModelsBlock)

    from Blocks.SetupGlide import setupGlideBlock  # type: ignore

    eapmPlugin.addBlock(setupGlideBlock)

    from Blocks.PeleEAPM import peleBlock  # type: ignore

    eapmPlugin.addBlock(peleBlock)

    from Blocks.AnalysePELE import analysePELEBlock  # type: ignore

    eapmPlugin.addBlock(analysePELEBlock)

    from Blocks.ConservedResiduesMSA import conservedResiduesMSABlock  # type: ignore

    eapmPlugin.addBlock(conservedResiduesMSABlock)

    from Blocks.Mafft import multipleSequenceAlignmentBlock  # type: ignore

    eapmPlugin.addBlock(multipleSequenceAlignmentBlock)

    from Blocks.HmmBuild import hmmBuildBlock  # type: ignore

    eapmPlugin.addBlock(hmmBuildBlock)

    from Blocks.HmmSearch import hmmSearchBlock  # type: ignore

    eapmPlugin.addBlock(hmmSearchBlock)

    from Blocks.HmmScan import hmmScanBlock  # type: ignore

    eapmPlugin.addBlock(hmmScanBlock)

    from Blocks.HmmAlign import hmmAlignBlock  # type: ignore

    eapmPlugin.addBlock(hmmAlignBlock)

    from Blocks.JackHmmer import jackHmmerBlock  # type: ignore

    eapmPlugin.addBlock(jackHmmerBlock)

    from Blocks.AsiteDesign import asiteDesignBlock  # type: ignore

    eapmPlugin.addBlock(asiteDesignBlock)

    # Add the configs
    from Configs.mafftConfig import mafftExecutableConfig  # type: ignore

    eapmPlugin.addConfig(mafftExecutableConfig)

    from Configs.hmmerConfig import hmmerExecutableConfig  # type: ignore

    eapmPlugin.addConfig(hmmerExecutableConfig)

    # Return the plugin
    return eapmPlugin


plugin = createPlugin()
