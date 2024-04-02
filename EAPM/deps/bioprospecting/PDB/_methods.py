from Bio import PDB
import os
import shutil

def retrievePDBs(pdb_codes, names=None, pdb_directory='PDB'):
    """
    Download a set of pdb structures from the PDB database.

    Parameters
    ----------
    pdb_codes : list or str
        A list with the PDB codes to retrieve their files. Optionally a string
        with a unique code can be given.
    names : dict
        Dictionary mapping the name of each PDB structure
    pdb_directory : str
        The name of the directory to store the retrieved files.

    Returns
    -------
    pdb_paths : list
        A list containing the paths to the retrieved PDB files.
    """

    pdb_paths = {}
   # get_alignment_length
    pdbl = PDB.PDBList()

    # Create directory to store files
    if not os.path.exists(pdb_directory):
        os.mkdir(pdb_directory)

    # Convert to list if only a string is given
    if isinstance(pdb_codes, str):
        pdb_codes = [pdb_codes]

    if isinstance(pdb_codes, list):
        # Iterate pdb codes
        for code in pdb_codes:
            # Download PDB file
            if names == None:
                output_file = pdb_directory+'/'+code.upper()+'.pdb'
            else:
                output_file = pdb_directory+'/'+names[code.upper()]+'.pdb'

            if not os.path.exists(output_file):
                pdbl.retrieve_pdb_file(code, file_format='pdb', pdir=pdb_directory)
            else: # If file already exists
                print('Structure exists: '+code.upper()+'.pdb')
                if names != None:
                    print('It was named as: '+names[code.upper()]+'.pdb')
                pdb_paths[code] = output_file

    for f in os.listdir(pdb_directory):
        if f.endswith('.ent'):
            # Rename file
            if names == None:
                renamed_file = pdb_directory+'/'+f.replace('pdb','').upper().replace('.ENT','.pdb')
            else:
                renamed_file = pdb_directory+'/'+names[f.replace('pdb','').upper().replace('.ENT','')]+'.pdb'
            os.rename(pdb_directory+'/'+f, renamed_file)
            # Append path
            pdb_paths[code] = renamed_file

    # Remove unnecesary folders created by Bio.PDB method
    shutil.rmtree('obsolete')

    return pdb_paths

def readPDB(pdb_file):
    """
    Reads a PDB file with the Biopython PDB parser.

    Parameters
    ----------
    pdb_file : str
        path to the pdb file.

    Returns
    -------
    structure : list or Bio.PDB.Structure
        Structure objects
    """

    parser = PDB.PDBParser()
    name = pdb_file.split('/')[-1].split('.pdb')[0]
    structure = parser.get_structure(name, pdb_file)

    return structure


def saveStructureToPDB(structure, output, remove_hydrogens=False, remove_water=False,
                        only_protein=False):
    """
    Saves a structure into a PDB file

    Parameters
    ----------
    structure : list or Bio.PDB.Structure
        Structure to save
    """

    io = PDB.PDBIO()
    io.set_structure(structure)

    selector = None
    if remove_hydrogens:
        selector = notHydrogen()
    elif remove_water:
        selector = notWater()
    elif only_protein:
        selector = onlyProtein()

    if selector != None:
        io.save(output, selector)
    else:
        io.save(output)

def getChainSequence(chain):
    """
    Get the one-letter protein sequence of a Bio.PDB.Chain object.

    Parameters
    ----------
    chain : Bio.PDB.Chain
        Input chain to retrieve its sequence from.

    Returns
    -------
    sequence : str
        Sequence of the input protein chain.
    None
        If chain does not contain protein residues.
    """
    sequence = ''
    for r in chain:
        if r.id[0] == ' ': # Non heteroatom filter
            try:
                sequence += PDB.Polypeptide.three_to_one(r.resname)
            except:
                sequence += 'X'
    if sequence == '':
        return None
    else:
        return sequence

def chainsAsStructure(chains):
    """
    This method creates a new Structure object containing only the given chains.

    Parameters
    ----------
    chains : list or Bio.PDB.Chain
        Chain or chains to be added to the new structure object.

    Returns
    -------
    structure : Bio.PDB.Structure
    """

    if not isinstance(chains, list):
        chains = [chains]

    structure = PDB.Structure.Structure(0)
    model = PDB.Model.Model(0)
    for chain in chains:
        model.add(chain)
    structure.add(model)

    return structure

def renumberResidues(structure, by_chain=False):
    """
    Renumber residues in a structure object starting from one. Two methods are possible:
    if by_chain is set to True the renumbering is restarted at the begining of every
    chain.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Input structure object
    by_chain : bool
        Whether each chain should be renumerated from one.

    Returns
    -------

    structure_copy : Bio.PDB.Structure
        Renumbered copy of the input structure
    """

    count = 0
    structure_copy = PDB.Structure.Structure(0)
    model = PDB.Model.Model(0)
    auxiliar = PDB.Chain.Chain(0)

    for chain in structure.get_chains():
        new_chain = PDB.Chain.Chain(chain.id)
        if by_chain:
            count = 0
        for residue in chain.get_residues():
            count += 1
            residue.set_parent(auxiliar)
            residue.id = (residue.id[0], count, residue.id[2])
            new_chain.add(residue)
        model.add(new_chain)
    structure_copy.add(model)

    return structure_copy
