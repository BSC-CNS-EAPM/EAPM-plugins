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

    from Blocks.AnalyseGlideDocking import analyseGlideDocking

    eapmPlugin.addBlock(analyseGlideDocking)

    from Blocks.PDBToMAE import convertPDBToMAEBlock

    eapmPlugin.addBlock(convertPDBToMAEBlock)

    from Blocks.TrimAlphafoldModels import trimAlphaFoldModelsBlock

    eapmPlugin.addBlock(trimAlphaFoldModelsBlock)

    from Blocks.SetupGlide import setupGlideBlock

    eapmPlugin.addBlock(setupGlideBlock)

    from Blocks.PeleEAPM import peleBlock

    eapmPlugin.addBlock(peleBlock)

    from Blocks.AnalysePELE import analysePELEBlock

    eapmPlugin.addBlock(analysePELEBlock)

    from Blocks.ConservedResiduesMSA import conservedResiduesMSABlock

    eapmPlugin.addBlock(conservedResiduesMSABlock)
    
    from Blocks.Mafft import multipleSequenceAlignmentBlock

    eapmPlugin.addBlock(multipleSequenceAlignmentBlock)
    
    from Blocks.MSA2HMM import convertMSAToHMMBlock

    eapmPlugin.addBlock(convertMSAToHMMBlock)
    
    from Blocks.HmmSearch import hmmsearchBlock

    eapmPlugin.addBlock(hmmsearchBlock)
    
    from Blocks.AsiteDesign import asiteDesignBlock

    eapmPlugin.addBlock(asiteDesignBlock)

    # Add the configs
    from Configs.mafftConfig import mafftExecutableConfig

    eapmPlugin.addConfig(mafftExecutableConfig)

    # Return the plugin
    return eapmPlugin


plugin = createPlugin()
