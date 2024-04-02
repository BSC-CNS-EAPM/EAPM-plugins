from .. import alignment
import os
import pandas as pd

def setNetsolpCalculations(job_folder, sequences, netsolp_path=None, n_splits=1, model='ESM1b', prediction='SU',
                           conda_env=None, check_unfinished=False):

    if netsolp_path == None:
        netsolp_script = '_netsolp_path/predict.py'
        models_path = '_netsolp_path/models'
    else:
        netsolp_script = netsolp_path+'/predict.py'
        models_path = netsolp_path+'/models'

    models = ['ESM12', 'ESM1b', 'Distilled', 'Both']
    if model not in models:
        raise ValueError('The given model type is not recognised. It should be %s' % str(models))

    predicitons = ['S', 'U', 'SU']
    if prediction not in predicitons:
        raise ValueError('The given prediction type is not recognised. It should be %s' % str(predicitons))

    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Partition predictions into splits

    # Create split folder
    split_size = int(len(sequences)/n_splits)
    remainder = len(sequences)-(split_size*n_splits)
    split_folder = {}
    splits = {}
    split_limit = {}
    for s in range(1, n_splits+1):
        split_folder[s] = job_folder+'/'+str(s).zfill(len(str(n_splits)))
        if not os.path.exists(split_folder[s]):
            os.mkdir(split_folder[s])
        splits[s] = {}
        split_limit[s] = split_size
        if s <= remainder:
            split_limit[s] += 1

    # Put sequences in different splits
    current = 1
    for s in sequences:
        splits[current][s] = sequences[s]
        if len(splits[current]) == split_limit[current]:
                alignment.writeFastaFile(splits[current],
                                         split_folder[current]+'/sequences.fasta')
                current += 1

    # Create execution script
    with open(job_folder+'/predict.sh', 'w') as nsps:
        if conda_env != None:
            nsps.write('eval "$(conda shell.bash hook)"\n')
            nsps.write('conda activate netsolp\n')
        nsps.write('python ')
        nsps.write(netsolp_script+' ')
        nsps.write('--FASTA_PATH sequences.fasta ')
        nsps.write('--OUTPUT_PATH prediction.csv ')
        nsps.write('--MODEL_TYPE '+model+' ')
        nsps.write('--PREDICTION_TYPE '+prediction+' ')
        nsps.write('--MODELS_PATH '+models_path+' ')
        nsps.write('\n')
        if conda_env != None:
            nsps.write('conda deactivate\n')

    if check_unfinished:
        unfinished = readNetsolpResults(job_folder, return_unfinished=True, verbose=False)
    else:
        unfinished = split_folder.keys()

    jobs = []
    for s in unfinished:
        command = 'cd '+split_folder[s]+'\n'
        if netsolp_path == None:
            command += 'sed -i "s/_netsolp_path/NETSOLP_PATH/g" ../predict.sh\n'
        command += 'bash ../predict.sh\n'
        command += 'cd '+'../'*len(split_folder[s].split('/'))+'\n'
        jobs.append(command)

    return jobs

def readNetsolpResults(job_folder, return_unfinished=False, verbose=True, include_sequence=True):
    """
    Read NetsolP predicitons.

    Parameters
    ==========
    job_folder : str
        Path to the NetsolP calculation (see setNetsolpCalculations())
    return_unfinished : bool
        Check and return unfinished splits.
    include_sequence : bool
        Include a column with the entry sequence in the output dataframe?
    """
    predictions = {}
    predictions['ID'] = []
    if include_sequence:
        predictions['Sequence'] = []

    if return_unfinished:
        unfinished = []

    for split in sorted(os.listdir(job_folder)):

        # Check if it is a split folder
        try:
            int(split)
        except:
            continue

        # Check if split calculation has finished
        prediction_file = job_folder+'/'+split+'/prediction.csv'
        if not os.path.exists(job_folder+'/'+split+'/prediction.csv'):
            if verbose:
                print(f'Split {split} has not finished!')

            if return_unfinished:
                unfinished.append(int(split))

            continue

        if not return_unfinished:

            # Read predicition data
            data = pd.read_csv(prediction_file)

            # Store predicted values
            predictions['ID'] += data['sid'].tolist()
            if include_sequence:
                predictions['Sequence'] += data['fasta'].tolist()
            if 'predicted_solubility' in data:
                predictions.setdefault('Solubility',[])
                predictions['Solubility'] += data['predicted_solubility'].tolist()
            if 'predicted_usability' in data:
                predictions.setdefault('Usability',[])
                predictions['Usability'] += data['predicted_usability'].tolist()

    if return_unfinished:
        return unfinished

    # Convert to dataframe
    predictions = pd.DataFrame(predictions).set_index('ID')

    return predictions
