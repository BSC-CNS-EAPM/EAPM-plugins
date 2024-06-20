from .. import alignment
import pandas as pd
import subprocess
import os
import shutil

def calculateSoluprotScores(target_sequences, program=None, output_file=None, overwrite=False):
    """
    Get the soluprot values for a set of target sequences.

    Parameters
    ==========
    target_sequences : dict
        Target sequences with sequence IDs as keys and string sequences as values.
    program : str
        Path to the SoluProt executable
    """

    # Check if executable path was given.
    if program == None:
        raise ValueError('You need to specify the path to the soluprot.py executable with program=aboslute_path!')
    if output_file == None or not isinstance(output_file, str):
        raise ValueError('You need to specify the path to the output file as a string!')
    if not output_file.endswith('.csv'):
        output_file += '.csv'

    if not os.path.exists(output_file) or overwrite:
        # Write sequence files
        alignment.writeFastaFile(target_sequences, '.seqs.tmp.fasta')

        # Run soluprot prediction
        with open('.soluprot.sh', 'w') as  spf:
            spf.write('source activate soluprot\n')
            spf.write('python '+program+' ')
            spf.write('--i_fa .seqs.tmp.fasta ')
            spf.write('--o_csv '+output_file+' ')
            spf.write('--tmp_dir .soluprot_tmp\n')
            spf.write('conda deactivate\n')
        subprocess.run('bash .soluprot.sh', shell=True)

        # Read soluprot result
        soluprot_results = pd.read_csv(output_file)
        soluprot_results.pop('runtime_id')
        soluprot_results.columns = ['Sequence', 'Soluble']
        soluprot_results.set_index('Sequence', inplace=True)
        soluprot_results.to_csv(output_file)

        # Remove input and output
        shutil.rmtree('.soluprot_tmp')
        os.remove('.soluprot.sh')
        os.remove('.seqs.tmp.fasta')

    else:
        soluprot_results = pd.read_csv(output_file)
        soluprot_results.set_index('Sequence', inplace=True)

    return soluprot_results

def calculateSoluprotParallel(job_folder, all_sequences, program=None,
                              overwrite=False, batches=1000):
    """
    Get the soluprot values for a set of target sequences.

    Parameters
    ==========
    job_folder : str
        Folder where to store the calculation.
    target_sequences : dict
        Target sequences with sequence IDs as keys and string sequences as values.
    program : str
        Path to the SoluProt executable
    """

    # Check if executable path was given.
    if program == None:
        raise ValueError('You need to specify the path to the soluprot.py executable with program=aboslute_path!')

    all_target_sequences = all_sequences

    jobs = []

    # Create job folder
    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    seq_per_batch = len(all_target_sequences.keys())/batches

    for batch in range(1,batches+1):

        batch_num = str(batch)
        if not os.path.exists(job_folder+'/'+'sequences_'+batch_num+'/'):
            os.mkdir(job_folder+'/'+'sequences_'+batch_num+'/')
            target_sequences = {}
            used_seqs = []
            for i,j in enumerate(all_target_sequences):
                target_sequences[j] = all_target_sequences[j]
                used_seqs.append(j)
                if i > seq_per_batch-1:
                    break
            for seq in used_seqs:
                all_target_sequences.pop(seq)

            alignment.writeFastaFile(target_sequences, job_folder+'/'+'sequences_'+batch_num+'/sequences_'+batch_num+'.fasta')

        # Create commands
        command = 'source activate soluprot\n'
        command += 'cd '+job_folder+'/'+'sequences_'+batch_num+'/\n'
        command += 'python '+program+' '
        command += '--i_fa '+'sequences_'+batch_num+'.fasta'+' '
        command += '--o_csv '+'output'+batch_num+'.csv'+' '
        command += '--tmp_dir '+'tmp_'+batch_num+' '
        command += '--no_tmhmm'+'\n'
        command += 'cd ../../\n'
        command += 'source deactivate\n'
        jobs.append(command)

    return jobs
