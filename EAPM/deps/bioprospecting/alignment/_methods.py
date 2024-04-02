from Bio import SeqIO, AlignIO
import json
import numpy as np
from scipy.spatial import distance_matrix
from Bio.PDB.Polypeptide import three_to_one
from . import _mafft_functions
from . import _methods
from multiprocessing import Pool, cpu_count
import tqdm
import os

def readFastaFile(fasta_file):
    """
    Read a fasta file and get the sequences as a dictionary

    Parameters
    ----------
    fasta_file : str
        Path to the input fasta file

    Returns
    -------
    sequences : dict
        Dictionary containing the IDs and squences in the fasta file.
    """

    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    return sequences

def writeFastaFile(sequences, output_file):
    """
    Write sequences to a fasta file.

    Parameters
    ----------
    sequences : dict
        Dictionary containing as values the strings representing the sequences
        of the proteins to align and their identifiers as keys.

    output_file : str
        Path to the output fasta file
    """

    # Write fasta file containing the sequences
    with open(output_file, 'w') as of:
        for name in sequences:
            of.write('>'+name+'\n')
            of.write(sequences[name]+'\n')

def writeMsaToFastaFile(msa, output_file,exclude=None, plain_sequences=False):
    """
    Write sequences inside an MSA to a fasta file.

    Parameters
    ----------
    msa : Bio.AlignIO
        Multiple sequence aligment in Biopython format.

    output_file : str
        Path to the output fasta file

    exclude : list
        list of msa id to exclude from output file
    """

    if exclude == None:
        exclude = []
    # Write fasta file containing the sequences
    with open(output_file, 'w') as of:
        for s in msa:
            if s.id  not in exclude:
                of.write('>'+s.id+'\n')
                if plain_sequences:
                    of.write(str(s.seq).replace('-','')+'\n')
                else:
                    of.write(str(s.seq)+'\n')


def readMsaFromFastaFile(alignment_file):
    """
    Read an MSA from a fasta file.

    Parameters
    ----------
    alignment_file : str
        Path to the alignment fasta file

    msa : Bio.AlignIO
        Multiple sequence aligment in Biopython format.
    """

    msa = AlignIO.read(alignment_file, 'fasta')

    return msa

def savePSIBlastAsJson(psi_blast_result, output_file):
    """
    Save the results of a psiblast calculation (e.g., _blast_functions.PSIBlastDatabase())
    into a json file.

    Parameters
    ==========
    psi_blast_result : dict
        Output dictionary from PSIBlastDatabase() function.
    output_file : str
        Path to the output file.

    Returns
    =======
    output_file : str
        Path to the output file.
    """
    with open(output_file, 'w') as of:
        json.dump(psi_blast_result, of)
    return output_file

def readPSIBlastFromJson(json_file):
    """
    Save the results of a psiblast calculation (e.g., _blast_functions.PSIBlastDatabase())
    into a json file.

    Parameters
    ==========
    json_file : str
        Path to the json file contaning the blas result

    Returns
    =======
    psi_blast_result : dict
        Result from the PSI blast calculation
    """
    with open(json_file) as jf:
        psi_blast_result = json.load(jf)
    return psi_blast_result


def msaIndexesFromSequencePositions(msa, sequence_id, sequence_positions):
    """
    Get the multiple sequence alignment position indexes matching those positions (zero-based) of a specific target sequence.

    Parameters
    ==========
    msa : Bio.AlignIO
        Multiple sequence aligment in Biopython format.
    sequence_id : str
        ID of the target sequence
    sequence_positions : list
        Target sequence positions to match (one-based indexes)

    Returns
    =======
    msa_indexes : list
        MSA indexes matching the target sequence positions (zero-based indexes)
    """

    if isinstance(sequence_positions, int):
        sequence_positions = [sequence_positions]

    msa_indexes = []
    p = 0

    # Check whether the given ID is presetnin the MSA.
    msa_ids = [x.id for x in msa]
    if sequence_id not in msa_ids:
        raise ValueError('Entry %s not found in MSA' % sequence_id)

    # Gather MSA index positions mathing the target sequence positions.
    updated = False
    for i in range(msa.get_alignment_length()):
        for a in msa:
            if a.id == sequence_id:
                if a.seq[i] != '-':
                    p += 1
                    updated = True
                else:
                    updated = False
        if p in sequence_positions and updated:
            msa_indexes.append(i)

    return msa_indexes

def getSequencePositionsFromMSAindexes(msa, indexes, return_identity=False):
    """
    For the MSA indexes position return the sequence positions of each entry in the given MSA.
    If the entry has a '-' at that position, then return None. Optionally, the aminoacid
    identity can be returned with the option return_identity=True. Importantly,
    the sequence positions are sorted.

    Parameters
    ==========
    msa : Bio.AlignIO
        Multiple sequence aligment in Biopython format.
    index : int
        The MSA index
    return_identity : bool
        Return identities instead of sequence positions.
    """

    # Initialize a sequence counter per each sequence
    position = {}
    for entry in msa:
        position[entry.id] = 0

    sequence_positions = {}
    # Count the sequence positions
    for i in range(msa.get_alignment_length()):
        for entry in msa:
            sequence_positions.setdefault(entry.id, [])
            if entry.seq[i] != '-':
                position[entry.id] += 1

            if i in indexes:
                # Return None if gap is found at position
                if entry.seq[i] == '-':
                    sequence_positions[entry.id].append(None)
                elif return_identity:
                    # Return the amino acid identity of the corresponding position
                    sequence_positions[entry.id].append(entry.seq[i])
                else:
                    sequence_positions[entry.id].append(position[entry.id])

    return sequence_positions

def checkMotifConservation(reference_sequence, target_sequences, reference_indexes, output_file=None,
                           skip_failed=True, cpus=1, overwrite=False):
    """
    Align sequences individually to the target sequence to search for motif conservation.

    Parameters
    ==========
    reference_sequence : str
        Sequence containing the motif at the specified positions.
    target_sequences : dict
        Sequences to be compared against the reference_sequence for motif conservation
    reference_indexes : list
        Sequence indexes (one-based) where the motif appears in the reference_sequence
    output_file : str
        Optional path to a fasta file to store the sequences having the motif conserved.

    Returns
    =======
    conservation : dict
        Dictionary containing the boolean values if the motif is conserved or not in the given target_sequences.
    """

    if not all(reference_indexes[i] <= reference_indexes[i+1] for i in range(len(reference_indexes) - 1)):
        raise ValueError('reference_index should be an orderly list.')

    if output_file == None or (not os.path.exists(output_file) or overwrite):

        if cpus > cpu_count():
            print('Max CPU is %s, using them all' % cpu_count())
            cpus = cpu_count()
        pool = Pool(cpus)

        # Get the motif of each protein
        motifs = {}
        count = 0
        if cpus > 1:
            jobs = []
        for code in target_sequences:
            count += 1
            if cpus > 1:
                jobs.append((code, target_sequences[code], reference_sequence, reference_indexes, skip_failed))
            else:
                motifs[code] = _checkMotif(target_sequences[code], reference_sequence, reference_indexes, skip_failed)
                print('Processed %s of %s sequences' % (count, len(target_sequences)), end='\r')

        if cpus > 1:
            for code, motif in tqdm.tqdm(pool.imap_unordered(_checkMotifParalell, jobs), total=len(jobs)):
                motifs[code] = motif

        if output_file:
            with open(output_file, 'w') as mf:
                json.dump(motifs, mf)

    elif output_file != None and os.path.exists(output_file):
        with open(output_file) as mf:
            motifs = json.load(mf)

    return motifs

def _checkMotif(target_sequence, reference_sequence, indexes, skip_failed):

    seqs = {}
    seqs['reference'] = reference_sequence
    seqs['target'] = target_sequence

    if skip_failed:
        try:
            msa = _mafft_functions.mafft.multipleSequenceAlignment(seqs, stderr=False)
        except:
            return None
    else:
        msa = _mafft_functions.mafft.multipleSequenceAlignment(seqs, stderr=False)

    msa_indexes = msaIndexesFromSequencePositions(msa, 'reference', indexes)
    motif= getSequencePositionsFromMSAindexes(msa, msa_indexes, return_identity=True)
    positions = getSequencePositionsFromMSAindexes(msa, msa_indexes)

    assert positions['reference'] == indexes

    return motif['target']

def _checkMotifParalell(arguments):
    code, target_sequence, reference_sequence, indexes, skip_failed = arguments
    motif = _checkMotif(target_sequence, reference_sequence, indexes, skip_failed)
    return code, motif

def trimSequences(refence_sequence, sequences):
    """
    Cut sequence terminal regions that do not align with the reference sequence

    Parameters
    ==========
    refence_sequence : str
        Aminoacid reference sequence
    sequences : dict
        Dictionary containing the sequences to trim

    Returns
    =======
    trimmed_sequences : dict
        Dictionary containing the trimmed sequences

    """
    trimmed_sequences = {}
    for i,m in enumerate(sequences):

        print(f'Trimming sequence {i} of {len(sequences)}', end='\r')
        seqs = {}
        seqs['reference'] = refence_sequence
        seqs['target'] = sequences[m]
        msa = bioprospecting.alignment.mafft.multipleSequenceAlignment(seqs, stderr=False)
        for e in msa:
            if e.id == 'reference':
                # Get msa index fo first reference position
                for i in range(msa.get_alignment_length()):
                    if e.seq[i] != '-':
                        first = i
                        break
                # Get msa index fo first reference position
                for i in range(msa.get_alignment_length()-1,0, -1):
                    if e.seq[i] != '-':
                        last = i
                        break

        # Get target sequence trim indexes
        indexes = []
        for e in msa:
            if e.id == 'target':
                count = 0
                for i in range(msa.get_alignment_length()):
                    if i == first:
                        indexes.append(count)
                    if i == last:
                        indexes.append(count)
                    if e.seq[i] != '-':
                        count += 1

        assert len(indexes) == 2
        trimmed_sequences[m] = sequences[m][indexes[0]:indexes[1]]

    return trimmed_sequences

def getAllowedPositionsFromMSA(msa, target_positions, target_sequence, return_propensity=False):
    """
    Get the allowed AA for each target_sequence position in the MSA.

    Parameters
    ==========
    msa : Bio.Align.MultipleSeqAlignment
        Input multiple sequence alginment
    target_positions : list
        List with the index of the target sequence positions
    target_sequence : str
        MSA-ID of the target sequence.
    return_propensity : dict
        Return a nested dictionary containing the target positions and the allowed AAs as keys and
        the the AA propensity as values.
    Returns
    =======
    mutatable_positions_aas : dict
        Dictionary containing the target positions as keys and the allowed AAs at each as values.
    """

    count = {}

    allowed_aa = {}
    for i in range(msa.get_alignment_length()):
        for s in msa:
            # Get target sequence
            if s.id == target_sequence:
                t_sequence = s.seq
            count.setdefault(s.id, 0)
            if s.seq[i] != '-':
                count[s.id] += 1
            if count[target_sequence] in target_positions and t_sequence[i] != '-':
                position = count[target_sequence]
                allowed_aa.setdefault(position, {})
                allowed_aa[position].setdefault(s.seq[i], 0)
                allowed_aa[position][s.seq[i]] += 1

    # Remove non AA characters
    for p in allowed_aa:
        if '-' in allowed_aa[p]:
            allowed_aa[p].pop('-')
        if 'X' in allowed_aa[p]:
            allowed_aa[p].pop('X')
        if 'Z' in allowed_aa[p]:
            allowed_aa[p].pop('Z')

    if return_propensity:
        return allowed_aa

    # Compile dictionary
    mutatable_positions_aas = {}
    for p in allowed_aa:
        mutatable_positions_aas[p] = []
        for aa in allowed_aa[p]:
            mutatable_positions_aas[p].append(aa)
    return mutatable_positions_aas

def getAlignedResiduesBasedOnStructuralAlignment(ref_struct,target_struct,max_ca_ca=5.0,res_seq_gap=5):
    """
    Return a sequence string with aligned residues based on a structural alignment. All residues
    not structurally aligned are returned as '-'.
    Parameters
    ==========
    ref_struct : str
        Reference structure
    target_struct : str
        Target structure
    max_ca_ca : float
        Maximum CA-CA distance to be considered aligned
    Returns
    =======
    aligned_residues : str
        Full-length sequence string containing only aligned residues.
    """

    # Get sequences
    r_sequence = ''.join([three_to_one(r.resname) for r in ref_struct.get_residues() if r.id[0] == ' '])
    t_sequence = ''.join([three_to_one(r.resname) for r in target_struct.get_residues() if r.id[0] == ' '])

    msa = bioprospecting.alignment.mafft.multipleSequenceAlignment({'target':t_sequence,'reference':r_sequence})

    # Get alpha-carbon coordinates
    r_ca_coord = np.array([a.coord for a in ref_struct.get_atoms() if a.name == 'CA' and a.get_parent().id[0] == ' '])
    t_ca_coord = np.array([a.coord for a in target_struct.get_atoms() if a.name == 'CA' and a.get_parent().id[0] == ' '] )

    # Map residues based on a CA-CA distances
    D = distance_matrix(t_ca_coord, r_ca_coord) # Calculate CA-CA distance matrix
    D =  np.where(D <= max_ca_ca, D, np.inf) # Cap matrix to max_ca_ca

    # Map residues to the closest CA-CA distance
    mapping = {}

    # Start mapping from closest distance avoiding double assignments
    while not np.isinf(D).all():
        dists = np.unravel_index(np.argsort(D,axis=None), D.shape)
        dists = tuple(zip(dists[0],dists[1]))

        for i,j in dists:
            if D[i,j] != np.inf:
                r_msa_pos = bioprospecting.alignment.msaIndexesFromSequencePositions(msa,'reference',int(j)+1)
                t_msa_pos = bioprospecting.alignment.msaIndexesFromSequencePositions(msa,'target',int(i)+1)

                if abs(r_msa_pos[0] - t_msa_pos[0]) < res_seq_gap:
                    mapping[j] = i
                    D[i].fill(np.inf) # This avoids double assignments
                    D[:,j].fill(np.inf) # This avoids double assignments
                    break

                else:
                    D[i,j] = np.inf

                #else D[i,j] == np.inf:
                #    D[i].fill(np.inf) # This avoids double assignments
                #    D[:,j].fill(np.inf) # This avoids double assignments
                    #print('break')
                    #break
    # Create alignment based on the structural alignment mapping
    aligned_residues = []

    for i,r in enumerate(r_sequence):
        if i in mapping:
            #aligned_residues[model][(i+1,int(mapping[i])+1)] = (r,t_sequence[int(mapping[i])])
            aligned_residues.append(t_sequence[int(mapping[i])])
        else:
            #aligned_residues[model][(i+1,'-')] = (r,'-')
            aligned_residues.append('-')

    # Join list to get aligned sequences
    aligned_residues = ''.join(aligned_residues)

    return aligned_residues
