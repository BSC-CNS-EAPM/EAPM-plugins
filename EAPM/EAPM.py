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

    from Blocks.AlignPdbEAPM import alignBlock  # type: ignore

    eapmPlugin.addBlock(alignBlock)

    from Blocks.AnalyseAFConfidence import analyse_AF_confidence_block  # type: ignore

    eapmPlugin.addBlock(analyse_AF_confidence_block)
    
    from Blocks.AlignPDBMaestro import align_pdb_schro_block

    eapmPlugin.addBlock(align_pdb_schro_block)
    # Return the plugin
    return eapmPlugin


plugin = createPlugin()
