from .. import alignment
from .. import databases
from .pdb_data import pdb_list

from Bio import SeqIO
import io
import os
import numpy as np
from pkg_resources import resource_stream, Requirement

def getUniprotSequences(uniprot_codes):
    """
    Get the sequences for a list of Uniprot IDs.
    """
    if isinstance(uniprot_codes, str):
        uniprot_codes = [uniprot_codes]
    if not isinstance(uniprot_codes, list):
        raise ValueError('Variable uniprot_codes must be a string or a list!')

    uniprot_db = databases.downloadUniprot()

    sequences = {}
    for record in SeqIO.parse(uniprot_db, "fasta"):
        up_code = record.id.split('|')[1]
        if up_code in uniprot_codes:
            sequences[up_code] = str(record.seq)

    for code in uniprot_codes:
        if code not in sequences:
            print('Warning: code %s was not found in UniProt database.' % code)

    return sequences

def blastUniProtDatabase(target_sequence, max_target_seqs=500):

    uniprot_db = databases.downloadUniprot()

    blast_sequences = alignment.blast.blastDatabase(target_sequence, uniprot_db,
                                  max_target_seqs=max_target_seqs)

    return blast_sequences

def PSIBlastUniProtDatabase(target_sequence, num_iterations=5, max_target_seqs=500,
                            output_file=None):

    uniprot_db = databases.downloadUniprot()

    psiblast_sequences = alignment.blast.PSIBlastDatabase(target_sequence, uniprot_db,
                                                      max_target_seqs=max_target_seqs,
                                                      num_iterations=num_iterations,
                                                      output_file=output_file)

    return psiblast_sequences

def searchTriads(job_folder, pdb_path):

    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Get control script
    script_name = "search_catalytic_triad.py"
    control_script = resource_stream(Requirement.parse("bioprospecting"),
                                         "bioprospecting/scripts/"+script_name)
    control_script = io.TextIOWrapper(control_script)

    # Write control script to job folder
    with open(job_folder+'/._'+script_name, 'w') as sof:
        for l in control_script:
            sof.write(l)

    jobs = []
    for pdb in os.listdir(pdb_path):
        if pdb.endswith('.pdb'):
            command = 'cd '+job_folder+'\n'
            command += 'python ._'+script_name+' -p ../'+pdb_path+'/'+pdb+'\n'
            command += 'cd ..\n'
            jobs.append(command)
    return jobs

def chainDictionary(job_folder, pdb_path):

    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Get control script
    script_name = "pdb_chain_dictionary.py"
    control_script = resource_stream(Requirement.parse("bioprospecting"),
                                         "bioprospecting/scripts/"+script_name)
    control_script = io.TextIOWrapper(control_script)

    # Write control script to job folder
    with open(job_folder+'/._'+script_name, 'w') as sof:
        for l in control_script:
            sof.write(l)

    jobs = []
    for pdb in pdb_list:
        command = 'cd '+job_folder+'\n'
        command += 'python ._'+script_name+' -p '+pdb_path+'/'+pdb+'\n'
        command += 'cd ..\n'
        jobs.append(command)

    return jobs

def triadClassification(job_folder):

    if not os.path.exists(job_folder):
        os.mkdir(job_folder)

    # Get control script
    script_name = "classify_found_triads.py"
    control_script = resource_stream(Requirement.parse("bioprospecting"),
                                         "bioprospecting/scripts/"+script_name)
    control_script = io.TextIOWrapper(control_script)

    # Write control script to job folder
    with open(job_folder+'/._'+script_name, 'w') as sof:
        for l in control_script:
            sof.write(l)

    jobs = []

    command = 'python '+job_folder+'/._'+script_name+'\n'
    jobs.append(command)

    return jobs
