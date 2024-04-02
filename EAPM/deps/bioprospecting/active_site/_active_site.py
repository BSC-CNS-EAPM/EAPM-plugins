import os
import shutil
import io
import mdtraj
import uuid
from pkg_resources import resource_stream, Requirement

import mdtraj as md
import numpy as np
from Bio.Align import substitution_matrices

def definePocketResiduesWithSiteMap(job_folder, pdb_files, target_residue, distance_to_points=2.5,
                                    pH=7.0, epik_pH=False, samplewater=False, epik_pHt=False, ligand=False,
                                    remove_hydrogens=True, delwater_hbond_cutoff=False, prepwizard_opt=True,
                                    fill_loops=False, site_box=10, resolution='fine', reportsize=100,
                                    skip_missing_commands=False, return_missing=False):
    """
    Run a site map calcualtion to define the residues in a specific protein pocket.

    Parameters
    ==========
    job_folder : str
        Path to the folder were the job will be run.
    pdb_files : list
        Path to the input PDB files
    target_residue : dict
        Target residue for each pdb file as (chain_id, residue_index)
    """

    if not isinstance(pdb_files, list):
        raise ValuError('You must give a list of PDB files.')

    # Create calculation folders
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    if not os.path.exists(job_folder+'/input_models'):
        os.mkdir(job_folder+'/input_models')

    if not os.path.exists(job_folder+'/output_models'):
        os.mkdir(job_folder+'/output_models')

    # Copy Schrodinger control script
    # _copyScriptFile(job_folder, 'schrodinger_control.py')

    # Copy site map script
    _copyScriptFile(job_folder, 'prepareForSiteMap.py')

    jobs = []
    residues = {}
    missing = []
    for pdb_path in pdb_files:

        # Setup input pdb paths
        file_name = pdb_path.split('/')[-1]
        pdb_name = file_name.replace('.pdb', '')

        if prepwizard_opt:
            input_protein = job_folder+'/output_models/'+pdb_name+'/'+file_name
        else:
            input_protein = job_folder+'/input_models/'+file_name


        # Calculate list of pocket residues points
        if os.path.exists(job_folder+'/output_models/'+pdb_name+'/'+pdb_name+'_site_1_volpts.pdb'):
            print('SiteMap output for %s found. Calculating pocket residues.' % pdb_name)

            vpt_protein_pdb = '.'+str(uuid.uuid4())+'.pdb'
            with open(vpt_protein_pdb, 'w') as pdb_and_volume:
                with open(input_protein) as protein:
                    for l in protein:
                        if l.startswith('ATOM') or l.startswith('HETATM'):
                            pdb_and_volume.write(l)

                with open(job_folder+'/output_models/'+pdb_name+'/'+pdb_name+'_site_1_volpts.pdb') as volume:
                    for l in volume:
                        if l.startswith('ATOM') or l.startswith('HETATM'):
                            pdb_and_volume.write(l)

            traj = md.load(vpt_protein_pdb)
            protein = traj.topology.select('protein and not resname vpt')
            vpts = traj.topology.select('resname vpt')
            n = md.compute_neighbors(traj, distance_to_points/10, vpts, haystack_indices=protein)
            residues[pdb_name] = set([traj.topology.atom(i).residue.resSeq for i in n[0] if traj.topology.atom(i).is_sidechain])
            os.remove(vpt_protein_pdb)
            residues[pdb_name] = np.array(list(residues[pdb_name]))

        else:
            print('SiteMap output for %s not found. Returning SiteMap commands.' % pdb_name)

            if return_missing:
                missing.append(pdb_name)

            if skip_missing_commands:
                continue

            # Copy input PDB
            shutil.copyfile(pdb_path, job_folder+'/input_models/'+file_name)
            ## Set up prepwizard calculation ##

            # Create model output folder
            output_folder = job_folder+'/output_models/'+pdb_name
            if not os.path.exists(output_folder):
                os.mkdir(output_folder)

            command = 'cd '+output_folder+'\n'

            # Prepwizard optimization
            if prepwizard_opt:
                command += '"${SCHRODINGER}/utilities/prepwizard" '
                command += '../../input_models/'+pdb_name+'.pdb '
                command += file_name+' '
                command += '-fillsidechains '
                command += '-disulfides '
                if remove_hydrogens:
                    command += '-rehtreat '
                if epik_pH:
                    command += '-epik_pH '+str(pH)+' '
                if epik_pHt:
                    command += '-epik_pHt '+str(epik_pHt)+' '
                command += '-propka_pH '+str(pH)+' '
                command += '-f 2005 '
                command += '-rmsd 0.3 '
                if samplewater:
                    command += '-samplewater '
                if delwater_hbond_cutoff:
                    command += '-delwater_hbond_cutoff '+str(delwater_hbond_cutoff)+' '
                command += '-JOBNAME prepare_'+pdb_name+' '
                command += '-HOST localhost:1 '
                command += '-WAIT\n'

            ## Set up site Map calculation ##

            # Add command to convert PDB into sitemap format
            command += '"${SCHRODINGER}/run" ../../._prepareForSiteMap.py '
            if prepwizard_opt:
                command += pdb_name+'.pdb '
            else:
                command += '../../input_models/'+pdb_name+'.pdb '
            command += '. '
            if not ligand:
                command += '--protein_only'
            command += '\n'

            # Add site map calculation command
            command += '"${SCHRODINGER}/sitemap" '
            if not ligand:
                command += '-j '+pdb_name+' '
                command += '-prot '+pdb_name+'_protein.mae'+' '
                command += '-sitebox '+str(site_box)+' '
                command += '-resolution '+str(resolution)+' '
                command += '-reportsize '+str(reportsize)+' '
                command += '-keepvolpts yes '
                command += '-keeplogs yes '
                command += '-siteasl \"chain.name {'+str(target_residue[pdb_name][0])+'} '
                command += 'res.num {'+str(target_residue[pdb_name][1])+'}\" '
            else:
                input_ligand = pdb_name+'_ligand.mae'
                input_protein = pdb_name+'_protein.mae'

                with open(job_folder+'/output_models/'+pdb_name+'/'+pdb_name+'.in', 'w') as smi:
                    smi.write('PROTEIN ../../input_models/'+input_protein+'\n')
                    smi.write('LIGMAE ../../input_models/'+input_ligand+'\n')
                    smi.write('SITEBOX '+str(site_box)+'\n')
                    smi.write('RESOLUTION '+resolution+'\n')
                    smi.write('REPORTSIZE 100\n')
                    smi.write('KEEPVOLPTS yes\n')
                    smi.write('KEEPLOGS yes\n')

                # Add site map command
                command += pdb_name+'.in'+' '

            command += '-HOST localhost:1 '
            command += '-TMPLAUNCHDIR '
            command += '-WAIT\n'
            command += 'cd ../../..\n'
            jobs.append(command)

    if return_missing:
        return missing
    elif jobs != []:
        return jobs
    else:
        return residues

def _copyScriptFile(output_folder, script_name):
    """
    Copy a script file from the prepare_proteins package.

    Parameters
    ==========

    """
    # Get control script
    control_script = resource_stream(Requirement.parse("bioprospecting"),
                                     "bioprospecting/scripts/"+script_name)
    control_script = io.TextIOWrapper(control_script)

    # Write control script to output folder
    with open(output_folder+'/._'+script_name, 'w') as sof:
        for l in control_script:
            sof.write(l)

def getSequenceSimilarityDistanceMatrix(models, models_residues, models_sequences, mode='union',
                                        substmatrix='BLOSUM62', gap_penalty=-5):
    """
    Calculate a sequence similarity distance matrix between the given models for a subset of residue.

    The values are obtained summing the BLOSUM scores for each aligned pair of residues in the subset
    of residues to compare. The scores are first divided by the length of the compared residues. This
    makes the scores independent of the number of residues compared in the active site. Then, the full
    matrix is normalized to lay between 0 and 1. Finally, values are converted to a distance by substracting
    them to 1.

    Parameters
    ==========
    models : list
        List of models to compare.
    models_residues : dict
        Dictionary containing the residues subsets to compare for each model.
    models_sequences : dict
        Dictionary containing the sequences of the models to compare.
    mode : str
        Compare all residues from both sequences (union) or only the common residues (intersection).
    substmatrix : str
        Name of the substitution matrix to employ to calculate similarity scores
        See Bio.Align.substitution_matrices for details.

    Returns
    =======
    M : np.ndarray
        Normalized sequence similarity matrix.
    """

    def compareSequences(mi, mj, mode='union'):
        """
        Compare two sequences at the specified residues using a Blossum matrix.
        """

        seqs = {}
        seqs[mi] = models_sequences[mi]
        seqs[mj] = models_sequences[mj]
        msa = bioprospecting.alignment.mafft.multipleSequenceAlignment(seqs, stderr=False)

        idxi = bioprospecting.alignment.msaIndexesFromSequencePositions(msa, mi, models_residues[mi])
        idxj = bioprospecting.alignment.msaIndexesFromSequencePositions(msa, mj, models_residues[mj])

        # Define a mode for selecting residues to be compared
        if mode == 'union':
            idx = list(set(idxi+idxj))
        elif mode == 'intersection':
            idx = [x for x in idxi if x in idxj]

        # Calculate similarity score for the selected msa indexes
        s = 0
        seqi = []
        seqj = []
        for x in range(msa.get_alignment_length()):
            if x in idx:
                for e in msa:
                    if e.id == mi:
                        pi = e.seq[x]
                        if mi == mj:
                            pj = e.seq[x]
                    elif e.id == mj:
                        pj = e.seq[x]
                if pi != '-' and pj != '-':
                    ai = blosum.alphabet.index(pi)
                    aj = blosum.alphabet.index(pj)
                    s += blosum[ai][aj]
                else:
                    s += gap_penalty
                seqi.append(pi)
                seqj.append(pj)

        # Normalise by number of compared residues

        s = s/len(idx)
        return s

    # Define substitution matrix
    blosum = substitution_matrices.load(substmatrix)

    modes = ['union', 'intersection']

    if mode not in modes:
        raise ValueError('The selected mode %s is not implemented. Please use one of these: %s or %s' % (mode, *modes))

    # Check that models have residues and sequences defined in the corresponding dictionaries
    for m in models:
        if m not in models_residues:
            raise ValueError('Model %s not found int models_residues dictionary')
        if m not in models_sequences:
            raise ValueError('Model %s not found int models_sequences dictionary')

    # Create empty similarity matrix
    M = np.zeros((len(models), len(models)))

    # Compare models at the given residues
    for i,mi in enumerate(models):
        for j,mj in enumerate(models):
            if i >= j:
                M[i][j] = compareSequences(mi, mj)
                M[j][i] = M[i][j] # Fill the symmetric diagonal

    # Set minimum value to zero
    M = M-np.min(M)
    # Convert similarity into distance
    M = 1-(M/np.max(M))

    return M
