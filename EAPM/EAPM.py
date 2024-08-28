"""
Entry point for the EAPM plugin
"""

from HorusAPI import Plugin


def create_plugin():
    """
    Generates the EAPM plugin and returns the instance
    """
    # ========== Plugin Definition ========== #

    eapm_plugin = Plugin(id="EAPM")

    # ========== Blocks ========== #
    # pylint: disable=import-outside-toplevel
    from Blocks.alphafold import alphafoldBlock

    # Add the block to the plugin
    eapm_plugin.addBlock(alphafoldBlock)

    from Blocks.prepwizard import prepWizardBlock

    eapm_plugin.addBlock(prepWizardBlock)

    from Blocks.align_pdb import alignBlock

    eapm_plugin.addBlock(alignBlock)

    from Blocks.setup_docking_grid import setupDockingGrid

    eapm_plugin.addBlock(setupDockingGrid)

    # from Blocks.analyse_glide_docking import analyseGlideDocking

    # eapm_plugin.addBlock(analyseGlideDocking)

    from Blocks.pdb_2_mae import convertPDBToMAEBlock

    eapm_plugin.addBlock(convertPDBToMAEBlock)

    from Blocks.trim_alphafold_models import trimAlphaFoldModelsBlock

    eapm_plugin.addBlock(trimAlphaFoldModelsBlock)

    from Blocks.setup_glide import setupGlideBlock

    eapm_plugin.addBlock(setupGlideBlock)

    from Blocks.pele import peleBlock

    eapm_plugin.addBlock(peleBlock)

    from Blocks.analyse_PELE import analysePELEBlock

    eapm_plugin.addBlock(analysePELEBlock)

    from Blocks.conserved_residues import conservedResiduesMSABlock

    eapm_plugin.addBlock(conservedResiduesMSABlock)

    from Blocks.mafft import multipleSequenceAlignmentBlock

    eapm_plugin.addBlock(multipleSequenceAlignmentBlock)

    from Blocks.hmmbuild import hmmBuildBlock

    eapm_plugin.addBlock(hmmBuildBlock)

    from Blocks.hmmsearch import hmmSearchBlock

    eapm_plugin.addBlock(hmmSearchBlock)

    from Blocks.hmmscan import hmmScanBlock

    eapm_plugin.addBlock(hmmScanBlock)

    from Blocks.hmmalign import hmmAlignBlock

    eapm_plugin.addBlock(hmmAlignBlock)

    from Blocks.jackhmmer import jackHmmerBlock

    eapm_plugin.addBlock(jackHmmerBlock)

    from Blocks.asitedesign import asiteDesignBlock

    eapm_plugin.addBlock(asiteDesignBlock)

    from Blocks.ahatool import ahatoolBlock

    eapm_plugin.addBlock(ahatoolBlock)

    from Blocks.ep_pred import epPredBlock

    eapm_plugin.addBlock(epPredBlock)

    # from Blocks.testBlock import testBlock

    # eapm_plugin.addBlock(testBlock)

    from Blocks.analyse_glide import AnalyseGBlock

    eapm_plugin.addBlock(AnalyseGBlock)

    from Blocks.Rbcavity import rbCavityBlock

    eapm_plugin.addBlock(rbCavityBlock)

    from Blocks.Rbdock import rbDockBlock

    eapm_plugin.addBlock(rbDockBlock)

    from Blocks.RbParameterFile import rbParameterFile

    eapm_plugin.addBlock(rbParameterFile)

    # Add the configs
    from Configs.mafftConfig import mafftExecutableConfig

    eapm_plugin.addConfig(mafftExecutableConfig)

    from Configs.hmmerConfig import hmmerExecutableConfig

    eapm_plugin.addConfig(hmmerExecutableConfig)

    from Pages.load_tables import load_page

    eapm_plugin.addPage(load_page)

    # pylint: enable=import-outside-toplevel

    # Return the plugin
    return eapm_plugin


plugin = create_plugin()
