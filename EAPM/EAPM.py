"""
Entry point for the EAPM plugin
"""

from HorusAPI import Plugin


def create_plugin():
    """
    Generates the EAPM plugin and returns the instance
    """
    # ========== Plugin Definition ========== #

    eapm_plugin = Plugin()

    # ========== Blocks (grouped by category) ========== #
    # pylint: disable=import-outside-toplevel

    # ---------- Structure Prediction ----------
    from Blocks.alphafold import alphafoldBlock

    eapm_plugin.addBlock(alphafoldBlock)

    from Blocks.trim_alphafold_models import trimAlphaFoldModelsBlock

    eapm_plugin.addBlock(trimAlphaFoldModelsBlock)

    # ---------- Structure Preparation ----------
    from Blocks.prepwizard import prepWizardBlock

    eapm_plugin.addBlock(prepWizardBlock)

    from Blocks.align_pdb import alignBlock

    eapm_plugin.addBlock(alignBlock)

    from Blocks.pdb_2_mae import convertPDBToMAEBlock

    eapm_plugin.addBlock(convertPDBToMAEBlock)

    from Blocks.collect_selected_pdbs import collectSelectedPDBsBlock

    eapm_plugin.addBlock(collectSelectedPDBsBlock)

    # ---------- Docking ----------
    from Blocks.setup_docking_grid import setupDockingGrid

    eapm_plugin.addBlock(setupDockingGrid)

    from Blocks.setup_glide import setupGlideBlock

    eapm_plugin.addBlock(setupGlideBlock)

    from Blocks.analyse_glide import AnalyseGBlock

    eapm_plugin.addBlock(AnalyseGBlock)

    # from Blocks.analyse_glide_docking import analyseGlideDocking

    # eapm_plugin.addBlock(analyseGlideDocking)

    # from Blocks.Rbcavity import rbCavityBlock

    # eapm_plugin.addBlock(rbCavityBlock)

    # from Blocks.Rbdock import rbDockBlock

    # eapm_plugin.addBlock(rbDockBlock)

    # ---------- PELE ----------
    from Blocks.pele import peleBlock

    eapm_plugin.addBlock(peleBlock)

    from Blocks.analyse_PELE import analysePELEBlock

    eapm_plugin.addBlock(analysePELEBlock)

    # ---------- Rosetta ----------
    from Blocks.rosetta_relax import rosettaRelaxBlock

    eapm_plugin.addBlock(rosettaRelaxBlock)

    from Blocks.analyse_rosetta_relax import analyseRosettaRelaxBlock

    eapm_plugin.addBlock(analyseRosettaRelaxBlock)

    # ---------- BioEmu ----------
    from Blocks.bioemu import bioEmuBlock

    eapm_plugin.addBlock(bioEmuBlock)

    from Blocks.analyse_bioemu import analyseBioEmuBlock

    eapm_plugin.addBlock(analyseBioEmuBlock)

    # ---------- Model Ranking ----------
    from Blocks.pareto_selection import paretoSelectionBlock

    eapm_plugin.addBlock(paretoSelectionBlock)

    # ---------- Sequence Alignment ----------
    from Blocks.mafft import multipleSequenceAlignmentBlock

    eapm_plugin.addBlock(multipleSequenceAlignmentBlock)

    from Blocks.conserved_residues import conservedResiduesMSABlock

    eapm_plugin.addBlock(conservedResiduesMSABlock)

    from Blocks.phylogenetic_tree import phylogeneticTreeBlock

    eapm_plugin.addBlock(phylogeneticTreeBlock)

    # ---------- Clustering & Selection ----------
    from Blocks.mmseqs2 import mmseqsClusterBlock

    eapm_plugin.addBlock(mmseqsClusterBlock)

    from Blocks.mmseqs2_slurm import mmseqsClusterSlurmBlock

    eapm_plugin.addBlock(mmseqsClusterSlurmBlock)

    from Blocks.mmseqs_threshold_sweep import mmseqsThresholdSweepBlock

    eapm_plugin.addBlock(mmseqsThresholdSweepBlock)

    from Blocks.select_cluster_representatives import selectClusterRepresentativesBlock

    eapm_plugin.addBlock(selectClusterRepresentativesBlock)

    from Blocks.combine_selected_sequences import combineSelectedSequencesBlock

    eapm_plugin.addBlock(combineSelectedSequencesBlock)

    from Blocks.sequence_length_distribution import sequenceLengthDistributionBlock

    eapm_plugin.addBlock(sequenceLengthDistributionBlock)

    # ---------- HMMER ----------
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

    # ---------- Design & Prediction ----------
    from Blocks.asitedesign import asiteDesignBlock

    eapm_plugin.addBlock(asiteDesignBlock)

    from Blocks.ahatool import ahatoolBlock

    eapm_plugin.addBlock(ahatoolBlock)

    from Blocks.ep_pred import epPredBlock

    eapm_plugin.addBlock(epPredBlock)

    from Blocks.codon_transformer import codonTransformerBlock

    eapm_plugin.addBlock(codonTransformerBlock)

    # from Blocks.testBlock import testBlock

    # eapm_plugin.addBlock(testBlock)

    # ---------- Documentation ----------
    from Blocks.prot_prot_documentation import prot_docs

    eapm_plugin.addBlock(prot_docs)

    # ========== Configs ========== #
    from Configs.mafftConfig import mafftExecutableConfig

    eapm_plugin.addConfig(mafftExecutableConfig)

    from Configs.hmmerConfig import hmmerExecutableConfig

    eapm_plugin.addConfig(hmmerExecutableConfig)

    from Configs.mmseqsConfig import mmseqsExecutableConfig

    eapm_plugin.addConfig(mmseqsExecutableConfig)

    from Configs.pyrosettaConfig import pyrosettaExecutableConfig

    eapm_plugin.addConfig(pyrosettaExecutableConfig)

    # ========== Pages ========== #
    from Pages.load_tables import load_page

    eapm_plugin.addPage(load_page)

    # pylint: enable=import-outside-toplevel

    # Return the plugin
    return eapm_plugin


plugin = create_plugin()
