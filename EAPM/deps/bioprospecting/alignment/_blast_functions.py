import os
import numpy as np
from collections import OrderedDict
import subprocess
import uuid
import json
from ._methods import *
from ._mafft_functions import *

class blast:
    """
    Class to hold methods to work with blast executable.

    Methods
    -------
    calculatePIDs()
        Fast method to calculate the PID of a sequence against many.
    _getPIDsFromBlastpOutput()
        Method to parse the ouput of blast and retrieve the PIDs.
    """

    def calculatePIDs(target_sequence, comparison_sequences):
        """
        Calculate the percentage of identity (PID) of a target sequence against a group
        of other sequences.

        Parameters
        ----------
        target_sequence : str
            Target sequence.
        comparison_sequences : list of str
            List of sequences to compare the target sequence against.

        Returns
        -------
        pids : numpy.array
            Array containg the PIDs of the target_sequence to all the comparison_sequences.
        """

        if isinstance(comparison_sequences, str):
            comparison_sequences = [comparison_sequences]
        elif not isinstance(comparison_sequences, list):
            raise ValueError('Comparison sequences must be given a single string or \
            as a list of strings.')

        # Write sequences into a temporary file
        with open('seq1.fasta.tmp', 'w') as sf1:
            sf1.write('>seq1\n')
            sf1.write(str(target_sequence))
        with open('seq2.fasta.tmp', 'w') as sf2:
            for i,s in enumerate(comparison_sequences):
                sf2.write('>seq'+str(i)+'\n')
                sf2.write(str(s)+'\n')
        # Execute blastp calculation
        try:
            os.system('blastp -query seq1.fasta.tmp -subject seq2.fasta.tmp -out ssblast.out.tmp -max_target_seqs '+str(len(comparison_sequences)))
        except:
            raise ValueError('blastp executable failed!')

        # Parse blastp output to extract PIDs
        pids = blast._getPIDsFromBlastpOutput('ssblast.out.tmp', len(comparison_sequences))
        pids = np.array(list(pids.values()))

        # Remove temporary files
        os.remove('seq1.fasta.tmp')
        os.remove('ssblast.out.tmp')
        os.remove('seq2.fasta.tmp')

        return pids

    def blastDatabase(job_folder, target_sequences, path_to_database_fasta, max_target_seqs=500):
        """
        Blast a specific sequence against a fasta database. The database fasta file
        must be supplied.
        """

        if not isinstance(target_sequences, dict):
            raise ValueError('target_sequences must be a dictionary of sequences!')

        if not os.path.exists(job_folder):
            os.mkdir(job_folder)

        max_target_seqs = str(max_target_seqs)

        blast_results = {}
        for ts in target_sequences:
            if not os.path.exists(job_folder+'/'+ts):
                os.mkdir(job_folder+'/'+ts)

            # Write target sequence into a temporary file
            input_fasta = job_folder+'/'+ts+'/'+ts+'.fasta'
            output_file = job_folder+'/'+ts+'/'+ts+'.out'

            if not os.path.exists(output_file):

                writeFastaFile({ts: target_sequences[ts]}, input_fasta)

                # Execute blastp calculation
                # Run blastp
                command = 'blastp '
                command += '-query '+input_fasta+' '
                command += '-subject '+path_to_database_fasta+' '
                command += '-out '+output_file+' '
                command += '-max_target_seqs '+max_target_seqs+' '
                try:
                    output = subprocess.check_output(command,stderr=subprocess.STDOUT,
                                                     shell=True, universal_newlines=True)
                except subprocess.CalledProcessError as exc:
                    print(exc.output, exc.returncode)
                    raise Exception("Problem with blastp execution.")

            blast_results[ts] = blast._parseBlastpOutput(output_file)

        return blast_results

    def PSIBlastDatabase(target_sequence, path_to_database_fasta, output_file=None, num_iterations=5,
                         max_target_seqs=500, output_pssm=False, overwrite=False):
        """
        Run a PSI BLAST search for a specific sequence against a fasta database.
        The database fasta file must be supplied.
        """

        if output_file != None and os.path.exists(output_file) and not overwrite:
            print('PSIBlast output file found. Reading results from %s' % output_file)
            with open(output_file) as of:
                blast_results = json.load(of)
        else:
            max_target_seqs = str(max_target_seqs)

            # Write target sequence into a temporary file
            target_file = '.'+str(uuid.uuid4())+'.fasta.tmp'
            tmp_file = '.'+str(uuid.uuid4())+'.out.tmp'
            with open(target_file, 'w') as sf1:
                sf1.write('>seq1\n')
                sf1.write(str(target_sequence))

            # Execute blastp calculation
            # Run psiblast
            command = 'psiblast -query '+target_file+' '
            command += '-subject '+path_to_database_fasta+' '
            command += '-out '+tmp_file+' '
            if output_pssm:
                command += '-out_pssm testpssm.smp  -save_pssm_after_last_round '
            command += '-max_target_seqs '+max_target_seqs+' '
            command += '-num_iterations '+str(num_iterations)+' '
            try:
                output = subprocess.check_output(command,stderr=subprocess.STDOUT,
                                                 shell=True, universal_newlines=True)

            except subprocess.CalledProcessError as exc:
                print(exc.output, exc.returncode)
                # Remove temporary files
                os.remove(target_file)
                os.remove(tmp_file)
                print("Problem with psiblast execution.")
                return None

            blast_results = blast._parsePSIBlastOutput(tmp_file)

            # Write results to output file if given
            if output_file != None:
                with open(output_file, 'w') as of:
                    json.dump(blast_results, of)

            # Remove temporary files
            os.remove(target_file)
            os.remove(tmp_file)

        return blast_results

    def createPSSM(sequences, output_file, db_file=None, stderr=True, stdout=True, save_msa=None):
        """
        Creates a PSSM based on the given sequences.

        Parameters
        sequences : (str, dict)
            Path to fasta file or a dictionary containing the input sequences.
        output_file : str
            Path to the output PSSM file.
        db_file : str
            Path to a blast database file for searching with the PSSM (output processing not implemented yet!)
        """
        # Manage stdout and stderr
        if stdout:
            stdout = None
        else:
            stdout = subprocess.DEVNULL

        if stderr:
            stderr = None
        else:
            stderr = subprocess.DEVNULL

        if isinstance(sequences, str):
            sequences = readFastaFile(sequences)

        first_seq_len = len(next(iter(sequences.values())))
        if not all(len(seq) == first_seq_len for seq in sequences.values()):
            # Create a multiple sequence alignment
            msa = mafft.multipleSequenceAlignment(sequences, stderr=False)
        else:
            print('All sequences are of the same length so no MSA calculation will be run.')
            msa = sequences

        if not save_msa:
            msa_file = '.'+str(uuid.uuid4())+'.fasta.tmp'
        else:
            msa_file = save_msa
        writeMsaToFastaFile(msa, msa_file)

        if not db_file:
            # Create a mock fasta file (necessary to run psiblast but calculation is not affected by it.)
            mock_file = '.'+str(uuid.uuid4())+'.fasta.tmp'
            with open(mock_file, 'w') as mf:
                mf.write('>mock_sequence\n')
                mf.write('AAA\n')
            db_file = mock_file

        # Create PSSM
        command =  'psiblast '
        command += '-db '+db_file+' '
        command += '-in_msa '+msa_file+' '
        command += '-out_ascii_pssm '+output_file+' '

        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)
        if not save_msa:
            os.remove(msa_file)
        os.remove(mock_file)

        return output_file

    def parsePSSM(pssm_file):
        """
        Parse PSSM file.

        Parameters
        ==========
        pssm_file : str
            Path to the PSSM file.

        Returns
        =======
        positions : numpy.ndarray
            PSSM positions of the first sequence in the MSA.
        aminoacids : numpy.ndarray
            The amino acids
        """
        log_llh = []
        positions = []
        native = []
        with open(pssm_file) as f:
            for l in f:
                # skip empty lines
                if l == '\n' or l.startswith('Last'):
                    continue
                # Get amino acids list
                elif l.split()[0] == 'A':
                    aminoacids = l.split()[:20]
                elif len(l.split()) == 44:
                    position, wt = l.split()[:2]
                    positions.append(int(position))
                    native.append(wt)
                    log_llh.append([int(x) for x in l.split()[2:22]])
                    l.split()[22:42]
                    l.split()[42:]

        aminoacids = np.array(aminoacids)
        positions = np.array(positions)
        log_llh = np.array(log_llh)

        return positions, native, aminoacids, log_llh

    def createLibraryFromPSSM(pssm_file, target_sequence=None, score_threshold=0, return_scores=False):
        """
        Create a mutational library based on the given PSSM log likelihood scores.
        """

        pssm = blast.parsePSSM(pssm_file)
        positions, native, aminoacids, log_llh = pssm
        # aminoacids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

        mutational_library = {}
        for p in range(1, log_llh.shape[0]+1):
            for i, aa in enumerate(aminoacids):
                if log_llh[p-1, i] >= score_threshold:
                    if return_scores:
                        mutational_library.setdefault(p, {})
                        mutational_library[p][aa] = int(log_llh[p-1, i])
                    else:
                        mutational_library.setdefault(p, [])
                        mutational_library[p].append(aa)

        return mutational_library

    def _parsePSIBlastOutput(psiblast_outfile):
        """
        Reads information from a psiblast output file.

        Parameters
        ----------
        psiblast_outfile : str
            Path to the psiblast outputfile.
        """
        sequences = {}

        # Analyse sequences by round
        # Read blast file and extract sequences full ids
        r = None
        with open(psiblast_outfile) as bo:
            cond = False
            for l in bo:
                if l.startswith('>'):
                    cond = True
                    full_name = l[2:][:-1]
                elif cond:
                    if l.startswith('Length='):
                        sequences[r].append(full_name)
                        cond = False
                    else:
                        full_name += l[:-1]
                if l.startswith('Results from round'):
                    r = int(l.split()[-1])
                    sequences[r] = []

        # Read blast file again to extract e-values matched to partial sequence names
        r = None
        blast_results = {}
        with open(psiblast_outfile) as bo:
            cond = False
            for l in bo:
                if 'Sequences producing significant alignments:' in l and r == 1:
                    cond = True
                elif 'Sequences used in model and found again:' in l:
                    continue
                elif 'Sequences not found previously or not previously below threshold:' in l:
                    if r > 1:
                        cond = True
                    continue
                elif l.startswith('>'):
                    cond = False
                elif cond and l.split() != []:
                    e_value = float(l.split()[-1])
                    name = l[:-25]
                    blast_results[r][name] = {}
                    blast_results[r][name]['e-value'] = e_value
                if l.startswith('Results from round'):
                    r = int(l.split()[-1])
                    blast_results[r] = {}

        # Match partial sequence names with full sequence names
        full_names = {}
        for round in blast_results:
            full_names[round] = {}
            for k1 in blast_results[round]:
                for k2 in sequences[round]:
                    if k1 in k2:
                        full_names[round][k1] = k2

        # Replace dict entries with full names
        for round in full_names:
            for k in full_names[round]:
                blast_results[round][full_names[round][k]] = blast_results[round].pop(k)

        return blast_results

    def _parseBlastpOutput(blastp_outfile, return_identity=True):
        """
        Reads information from blast output file.

        Parameters
        ----------
        blastp_outfile : str
            Path to the blastp outputfile.
        """
        sequences = []
        # Read blast file and extract sequences full ids
        with open(blastp_outfile) as bo:
            cond = False
            for l in bo:
                if l.startswith('>'):
                    cond = True
                    full_name = l[2:][:-1]
                elif cond:
                    if l.startswith('Length='):
                        sequences.append(full_name)
                        cond = False
                    else:
                        full_name += l[:-1]

        # Read blast file again to extract e-values matched to partial sequence names
        blast_results = {}
        with open(blastp_outfile) as bo:
            cond = False
            for l in bo:
                if 'Sequences producing significant alignments:' in l:
                    cond = True
                elif l.startswith('>'):
                    cond = False
                elif cond and l.split() != []:
                    e_value = float(l.split()[-1])
                    # name = l[:-20].strip()
                    name = l.split()[0]

                    blast_results[name] = {}
                    blast_results[name]['e-value'] = e_value

                if return_identity:
                    if l.startswith('>'):
                        code = l.split()[0].replace('>', '')
                        if code == '':
                            code = l.split()[1]
                    elif 'Identities' in l:
                        pid = float(l.split('(')[1].split(')')[0][:-1])/100
                        blast_results[code]['pid'] = pid

        # Match partial sequence names with full sequence names
        full_names = {}
        for k1 in blast_results:
            for k2 in sequences:
                if k1 in k2:
                    full_names[k1] = k2

        # Replace dict entries with full names
        for k in full_names:
            blast_results[full_names[k]] = blast_results.pop(k)

        return blast_results

    def _getPIDsFromBlastpOutput(blastp_outfile, n_sequences):
        """
        Parse output file from a blastp comparison to extract pids

        Parameters
        ----------
        blastp_outfile : str
            Path to the blastp outputfile.
        n_sequences : str
            number of sequences in the comparison file.

        Returns
        -------
        values : OrderedDict
            Dictionary containing the PID values.
        """

        # Create dictionary integer entries
        values = OrderedDict()
        for i in range(n_sequences):
            values[i] = 0

        # Read PID from blastp output file
        with open(blastp_outfile) as bf:
            for l in bf:
                if l.startswith('>'):
                    seq = int(l.split()[1].replace('seq',''))
                elif 'Identities' in l:
                    pid = eval(l.split()[2])
                    values[seq] = pid

        return values
