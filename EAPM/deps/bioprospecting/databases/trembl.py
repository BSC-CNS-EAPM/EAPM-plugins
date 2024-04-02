from .. import alignment
import os
import shutil
import Bio
import pandas as pd

def BLAST(job_folder, target_sequences, num_descriptions=5000, evalue=0.1, outfmt=None,
         max_target_seqs=None, get_sequences=True,
         trembl_path='/gpfs/projects/bsc72/uniprot/trembl/trembl'):
    """
    Do a BLAST search in the TREMBL database at a BSC cluster.

    Parameters
    ==========
    target_sequences : dict
        The sequences for which to run BLAST against TREMBL.
    num_descriptions : int
        The maximum number of outputs to get in the BLAST output.
    """
    if not isinstance(target_sequences, dict):
        raise ValueError('target_sequences must be a dictionary!')

    # Creat job folders
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    if get_sequences:
        if outfmt and outfmt != 6:
            print('Warning: when extrating sequences the outfmt is set to 6')
        outfmt = 6

    if num_descriptions != None and max_target_seqs != None:
        raise ValueError('Using num_descriptions and max_target_seqs is not compatible. Choose one of the two.')

    # Iterate each given sequence
    jobs = []
    for sequence in target_sequences:

        # Create sequence folder
        sdir = job_folder+'/'+sequence
        if not os.path.exists(sdir):
            os.mkdir(sdir)

        # Save fasta files in inputs folder
        alignment.writeFastaFile({sequence:target_sequences[sequence]},
                                 sdir+'/'+sequence+'.fasta')

        # Create command line for running PSI BLAST
        command = 'cd '+sdir+'\n'
        command += 'blastp '
        command += '-query '+sequence+'.fasta '
        command += '-db '+trembl_path+' '
        if num_descriptions != None:
            command += '-num_descriptions '+str(num_descriptions)+' '
        if max_target_seqs != None:
            command += '-max_target_seqs '+str(max_target_seqs)+' '
        if outfmt:
            command += '-outfmt '+str(outfmt)+' '
        command += '-out '+sequence+'.out '
        command += '-evalue '+str(evalue)+' '
        command += '\n'

        if get_sequences:
            command += "awk '{print $2}' "+sequence+".out > "+sequence+"_IDs.out\n"
            command += 'blastdbcmd '
            command += '-db '+trembl_path+' '
            command += '-entry_batch '+sequence+'_IDs.out '
            command += '> '+sequence+'.out.fasta\n'

        command += 'cd ../../\n'

        jobs.append(command)

    return jobs

def PSIBlast(job_folder, target_sequences, num_descriptions=5000, num_iterations=10, evalue=0.1,
             max_target_seqs=None, output_pssm=False, check_finished=True, outfmt=None,
             trembl_path = '/gpfs/projects/bsc72/uniprot/trembl/trembl', get_sequences=True,
             only_get_sequences=False, num_threads=None):
    """
    Do a PSIBLAST search in the TREMBL database at a BSC cluster.

    Parameters
    ==========
    target_sequences : dict

    """

    if get_sequences:
        if outfmt and outfmt != 6:
            print('Warning: when extrating sequences the outfmt is set to 6')
        outfmt = 6

    if outfmt in [6, 7, 10]:
        num_descriptions = None

    if not isinstance(target_sequences, dict):
        raise ValueError('target_sequences must be a dictionary!')

    # Creat job folders
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Iterate each given sequence
    jobs = []
    for sequence in target_sequences:

        # Create sequence folder
        sdir = job_folder+'/'+sequence
        if not os.path.exists(sdir):
            os.mkdir(sdir)

        if os.path.exists(sdir+'/'+sequence+'.out') and check_finished:
            blast_results = alignment.blast._parsePSIBlastOutput(sdir+'/'+sequence+'.out')
            if blast_results!= {} and max(blast_results.keys()) >= num_iterations:
                print('PSIBLAST for sequence %s found with %s iterations. Skipping calculation.' % (sequence, max(blast_results.keys())))
                continue

        # Save fasta files in inputs folder
        alignment.writeFastaFile({sequence:target_sequences[sequence]},
                                 sdir+'/'+sequence+'.fasta')

        if num_descriptions != None and max_target_seqs != None:
            raise ValueError('Using num_descriptions and max_target_seqs is not compatible. Choose one of the two.')

        # Create command line for running PSI BLAST

        command = 'cd '+sdir+'\n'

        if not only_get_sequences:
            command += 'psiblast '
            command += '-query '+sequence+'.fasta '
            if num_descriptions != None:
                command += '-num_descriptions '+str(num_descriptions)+' '
            if max_target_seqs != None:
                command += '-max_target_seqs '+str(max_target_seqs)+' '
            command += '-db '+trembl_path+' '
            command += '-out '+sequence+'.out '
            if output_pssm:
                command += '-out_ascii_pssm pssm.smp -save_each_pssm '
            if outfmt:
                command += '-outfmt '+str(outfmt)+' '
            command += '-num_iterations '+str(num_iterations)+' '
            command += '-evalue '+str(evalue)+' '
            if num_threads:
                command += '-num_threads '+str(num_threads)+' '
            command += '\n'

        if get_sequences or only_get_sequences:
            command += "awk '{print $2}' "+sequence+".out > "+sequence+"_IDs.out\n"
            command += 'blastdbcmd '
            command += '-db '+trembl_path+' '
            command += '-entry_batch '+sequence+'_IDs.out '
            command += '> '+sequence+'.out.fasta\n'

        command += 'cd '+'../'*len(sdir.replace('//','/').split('/'))+'\n'

        jobs.append(command)

    return jobs

def PSIBlastWithMSA(job_folder, msa, num_descriptions=5000, num_iterations=10, evalue=0.1,
                    max_target_seqs=None, output_pssm=False, check_finished=True, outfmt=None,
                    num_threads=1,
                    trembl_path = '/gpfs/projects/bsc72/uniprot/trembl/trembl', get_sequences=True,
                    only_get_sequences=False, negative_seqids=None):
    """
    Do a PSIBLAST search in the TREMBL database using multiple aligned sequences
    at a BSC cluster (or you can change the trembl_path keyword).

    Parameters
    ==========
    msa : Bio.Align.MultipleSeqAlignment

    negative_seqids : list
        List of seqids to exclude from the search.
    """

    if get_sequences:
        if outfmt and outfmt != 6:
            print('Warning: when extrating sequences the outfmt is set to 6')
        outfmt = 6

    if outfmt in [6, 7, 10]:
        num_descriptions = None

    if not isinstance(msa, Bio.Align.MultipleSeqAlignment):
        raise ValueError('msa must be a multiple sequence alignment as a Biopython Align MultipleSeqAlignment object!')

    # Creat job folders
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Create sequence folder
    sequence = 'MSA'

    sdir = job_folder+'/'+sequence
    if not os.path.exists(sdir):
        os.mkdir(sdir)

    if os.path.exists(sdir+'/'+sequence+'.out') and check_finished:
        blast_results = alignment.blast._parsePSIBlastOutput(sdir+'/'+sequence+'.out')
        if blast_results!= {} and max(blast_results.keys()) >= num_iterations:
            print('PSIBLAST for %s found with %s iterations. Skipping calculation.' % (sequence, max(blast_results.keys())))
            return

    # Save fasta files in inputs folder
    alignment.writeMsaToFastaFile(msa, sdir+'/'+sequence+'.fasta')

    if num_descriptions != None and max_target_seqs != None:
        raise ValueError('Using num_descriptions and max_target_seqs is not compatible. Choose one of the two.')

    if isinstance(negative_seqids, list):
        with open(sdir+'/negative_seqids.list', 'w') as lf:
            for code in negative_seqids:
                lf.write(code+'\n')
    elif negative_seqids:
        raise ValueError('negative_seqids should be a list!')

    # Create command line for running PSI BLAST
    command = 'cd '+sdir+'\n'

    if not only_get_sequences:
        command += 'psiblast '
        command += '-in_msa '+sequence+'.fasta '
        if num_descriptions != None:
            command += '-num_descriptions '+str(num_descriptions)+' '
        if max_target_seqs != None:
            command += '-max_target_seqs '+str(max_target_seqs)+' '
        command += '-db '+trembl_path+' '
        command += '-out '+sequence+'.out '
        if output_pssm:
            command += '-out_ascii_pssm pssm.smp -save_each_pssm '
        if outfmt:
            command += '-outfmt '+str(outfmt)+' '
        command += '-num_iterations '+str(num_iterations)+' '
        command += '-evalue '+str(evalue)+' '
        if num_threads > 1:
            command += '-num_threads '+str(num_threads)+' '
        if negative_seqids:
            command += '-negative_seqidlist negative_seqids.list '
        command += '\n'

    if get_sequences or only_get_sequences:
        command += "awk '{print $2}' "+sequence+".out > "+sequence+"_IDs.out\n"
        command += "sort "+sequence+"_IDs.out | uniq > "+sequence+"_IDs_uniq.out\n" # Remove redundancy if any
        command += "mv "+sequence+"_IDs_uniq.out "+sequence+"_IDs.out\n"
        command += 'blastdbcmd '
        command += '-db '+trembl_path+' '
        command += '-entry_batch '+sequence+'_IDs.out '
        command += '> '+sequence+'.out.fasta\n'

    command += 'cd '+'../'*len(sdir.replace('//','/').split('/'))+'\n'

    return [command]

def readBLASTResults(job_folder, only_codes=False):
    """
    Read the resuts from the PSIBlast() function over the TREMBL database.

    Parameters
    ==========
    job_folder : str
        Path to the job folder where the output of the PSIBlast() calculation is contained.
    only_codes : bool
        Split the sequence headers to only get the Uniprot codes.

    Returns
    =======
    sequences : dict
        The results from the PSIBLAST calculation.
    """

    sequences = {}
    for model in os.listdir(job_folder):

        model_out = job_folder+'/'+model+'/'+model+'.out'
        if not os.path.exists(model_out):
            print('BLAST output for model %s was not found!')
            continue
        else:
            blast_results = alignment.blast._parseBlastpOutput(model_out)
            if only_codes:
                sequences[model] = []
            else:
                sequences[model] = {}
            for s in blast_results:
                if only_codes:
                    sequences[model].append(s.split()[0])
                else:
                    sequences[model][s.split()[0]] = blast_results[s]['e-value']

    return sequences

def readPSIBLASTResults(job_folder, as_one_bundle=False, only_codes=False, remove_missing_folders=False,
                        verbose=True, outfmt=6):
    """
    Read the resuts from the PSIBlast() function over the TREMBL database.

    Parameters
    ==========
    job_folder : str
        Path to the job folder where the output of the PSIBlast() calculation is contained.
    as_one_bundle : bool
        Get the sequences as one unique list, i.e., without separation for iterations.
    only_codes : bool
        Split the sequence headers to only get the Uniprot codes.
    remove_missing_folders : bool
        Remove folders without output files.

    Returns
    =======
    sequences : dict
        The results from the PSIBLAST calculation.
    """

    implemented_outfmt = [6]

    if outfmt not in implemented_outfmt:
        raise ValueError('This output format has not been implemented yet.')

    if outfmt == 6:
        header = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

    frames = []
    for model in os.listdir(job_folder):

        model_out = job_folder+'/'+model+'/'+model+'.out'
        if not os.path.exists(model_out):
            if verbose:
                print('PSIBLAST output for model %s was not found!' % model)

            if remove_missing_folders:
                shutil.rmtree(job_folder+'/'+model)
            continue
        else:
            df = pd.read_csv(model_out,sep='\t',header=None)
            df.columns = header
            frames.append(df.set_index(['qseqid','sseqid'])[:-1])

    df = pd.concat(frames)

    return df

def readPSIBLASTSequences(job_folder, verbose=True):
    """
    Read the resuts from the PSIBlast() function over the TREMBL database.

    Parameters
    ==========
    job_folder : str
        Path to the job folder where the output of the PSIBlast() calculation is contained.

    Returns
    =======
    sequences : dict
        Sequences from the PSIBLAST calculation.
    """

    sequences = {}
    for model in os.listdir(job_folder):
        sequences[model] = {}
        fasta_out = job_folder+'/'+model+'/'+model+'.out.fasta'
        if not os.path.exists(fasta_out):
            if verbose:
                print('PSIBLAST output fasta was not found for model %s!' % model)

        else:
            model_sequences = alignment.readFastaFile(fasta_out)
            for code in model_sequences:
                if code not in sequences:
                    sequences[model][code] = model_sequences[code]

    return sequences
