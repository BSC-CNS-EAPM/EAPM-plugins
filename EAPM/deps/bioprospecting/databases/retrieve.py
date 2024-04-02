import os
import gzip
import shutil
import urllib
import uuid
import xml.etree.ElementTree as ET
from bioprospecting import alignment, databases, PDB
from multiprocessing import Pool, cpu_count

def getUniProtPDBChains(codes, output_folder, keep_original_pdbs=False):
    """
    Retrieve PDB chains associated with the specified uniprot codes.
    """

    # Create output folder to store original PDBs
    if not os.path.exists('PDB'):
        os.mkdir('PDB')

    # Create output folder to store PDB chains
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Get uniprot data to access PDB codes
    upd = databases.uniprot(codes)
    upd.getXMLFiles()
    upd_data = upd.parseXMLData()

    pdb_codes = {}
    for code in upd_data:

        # Create output folder to store original PDBs
        if not os.path.exists('PDB/'+code):
            os.mkdir('PDB/'+code)

        # Create output folder to store PDB chains
        if not os.path.exists(output_folder+'/'+code):
            os.mkdir(output_folder+'/'+code)

        # Download PDB structures
        pdb_codes[code] = upd_data[code]['Structures']

    # Scrape PDB database to get which chains are associated with the UP code.
    all_codes = []
    for code in pdb_codes:
        all_codes += pdb_codes[code]
    pdb_data = databases.scrape.getPDBInformation(all_codes, '.pdb.data')

    for code in upd_data:

        pdb_paths = PDB.retrievePDBs(pdb_codes[code], pdb_directory='PDB/'+code)

        # Store chains individually
        for pdb in pdb_paths:
            structure = PDB.readPDB(pdb_paths[pdb])
            chains = {c.id:c for c in structure.get_chains()}
            for c in pdb_data[pdb]['Macromolecules']:
                if pdb_data[pdb]['Macromolecules'][c]['Uniprot'] == code:
                    PDB.saveStructureToPDB(PDB.chainsAsStructure(chains[c]),
                                                          output_folder+'/'+code+'/'+pdb+'_'+c+'.pdb')

    # Remove folder with PDBS
    if not keep_original_pdbs:
        shutil.rmtree('PDB')

def createRepositoryFolder():
    """
    Create a hidden folder to store databases files.
    """
    home = os.path.expanduser("~")
    if not os.path.exists(home+'/.bioprospecting'):
        print('Library directory not found.')
        print('Creating library directory for the first time at %s' % home+'/.bioprospecting')
        os.mkdir(home+'/.bioprospecting')

    if not os.path.exists(home+'/.bioprospecting/databases'):
        os.mkdir(home+'/.bioprospecting/databases')

    database_dir = home+'/.bioprospecting/databases'

    return database_dir

def downloadUniprot(update=False):
    """
    Download Uniprot database
    """

    url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
    database_dir = createRepositoryFolder()

    if not os.path.exists(database_dir+'/uniprot'):
        os.mkdir(database_dir+'/uniprot')

    if not os.path.exists(database_dir+'/uniprot/uniprot_sprot.fasta') or update:
        if update:
            print('Updating UniProt database:')
        else:
            print('UniProt database not found.')
            print('Downloading UniProt database:')
        print('\tFrom '+url)
        print('\tAt '+database_dir+'/uniprot/uniprot_sprot.fasta')
        urllib.request.urlretrieve(url, database_dir+'/uniprot/uniprot_sprot.fasta.gz')

        f_in = gzip.open(database_dir+'/uniprot/uniprot_sprot.fasta.gz', 'rb')
        f_out = open(database_dir+'/uniprot/uniprot_sprot.fasta', 'wb')
        shutil.copyfileobj(f_in, f_out)
        f_in.close()
        f_out.close()
        os.remove(database_dir+'/uniprot/uniprot_sprot.fasta.gz')

    return database_dir+'/uniprot/uniprot_sprot.fasta'

def retrieveUniprotData(uniprot_code, database='uniprotkb'):
    """
    Retrieves UniProt data and stores it into a dictionary
    """

    available_dbs = ['uniprotkb', 'uniref']
    if database not in available_dbs:
        raise ValueError('Wrong database selected. Available databases %s ' % available_dbs)
    uniprot_data = {}

    # Retrieve sequence from fasta file
    if '_' in uniprot_code and database == 'uniprotkb':
        uniprot_code_search = uniprot_code.split('_')[0]
    else:
        uniprot_code_search = uniprot_code

    url = 'https://rest.uniprot.org/'+database+'/'+uniprot_code_search+'.fasta'
    fasta_file = '.'+str(uuid.uuid4())+'.fasta'

    try:
        urllib.request.urlretrieve(url, fasta_file)
    except:
        print('Problem downloading UniProt sequence for %s' % uniprot_code)
        print('Please check that the following URL exists %s' % url)

    if os.path.exists(fasta_file):
        sequence = alignment.readFastaFile(fasta_file)
        for code in sequence:
            uniprot_data['Sequence'] = sequence[code]
            break
    else:
        uniprot_data['Sequence'] = None

    if os.path.exists(fasta_file):
        os.remove(fasta_file)

    # Retrieve data from XML file

    url = 'https://rest.uniprot.org/'+database+'/'+uniprot_code_search+'.xml'
    xml_file = '.'+str(uuid.uuid4())+'.xml'

    try:
        urllib.request.urlretrieve(url, xml_file)
    except:

        print('Problem downloading UniProt information for %s' % uniprot_code)
        print('Please check that the following URL exists %s' % url)

    if not os.path.exists(xml_file):
        return None

    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Create dictionary entries
    uniprot_data['Protein'] = {}
    uniprot_data['Gene'] = {}
    uniprot_data['Organism'] = {}
    uniprot_data['GO'] = {}
    uniprot_data['Function'] = None
    uniprot_data['Reaction'] = []
    uniprot_data['Structures'] = []
    uniprot_data['features'] = {}

    tup = []

    # Iterate XML file to store UniProt data into a dictionary
    for x in root[0]:
        # Get sequence name
        if x.tag.endswith('name'):
            uniprot_data['Name'] = x.text

        # Get protein name information
        elif x.tag.endswith('protein'):
            for y in x:
                if y.tag.endswith('recommendedName'):
                    for z in y:
                        if 'fullName' == z.tag.split('}')[-1]:
                            uniprot_data['Protein']['recommendedName'] = z.text
                elif y.tag.endswith('alternativeName'):
                    for z in y:
                        if 'fullName' == z.tag.split('}')[-1]:
                            uniprot_data['Protein']['alternativeName'] = z.text

        # Get gene name information
        elif x.tag.endswith('gene'):
            for y in x:
                if y.tag.endswith('name'):
                    uniprot_data['Gene'][y.attrib['type']] = y.text

        # Get organism information
        elif x.tag.endswith('organism'):
            uniprot_data['Organism']['Name'] = {}
            for y in x:
                # Get organism names
                if 'type' in y.attrib:
                    uniprot_data['Organism']['Name'][y.attrib['type']] = y.text
                # Get taxonommic lineage
                elif 'lineage' == y.tag.split('}')[-1]:
                    uniprot_data['Organism']['Lineage'] = []
                    for taxon in y:
                        uniprot_data['Organism']['Lineage'].append(taxon.text)

        #Get features information
        elif x.tag.endswith('feature'):
            t = x.attrib['type']
            if t not in uniprot_data['features']:
                uniprot_data['features'][t] = []
            for y in x:
                if y.tag.endswith('location'):
                    for p in y:
                        if p.tag.endswith('position'):
                            pos = p.attrib['position']
                            uniprot_data['features'][t].append(pos)
                        elif p.tag.startswith("begin") or p.tag.startswith("end"):
                            pos = p.attrib['position']
                            if len(tup) == 1:
                                tup.append(pos)
                                uniprot_data['features'][t].append(tup)
                                tup = []
                            else:
                                tup.append(pos)

        # Process comments
        elif 'comment' == x.tag.split('}')[-1]:

            # Get function
            if 'type' in x.attrib and x.attrib['type'] == 'function':
                for y in x:
                    if 'text' == y.tag.split('}')[-1]:
                        uniprot_data['Function'] = y.text

            # Get catalytic activity
            if 'type' in x.attrib and x.attrib['type'] == 'catalytic activity':
                for y in x:
                    if y.tag.split('}')[-1] == 'reaction':
                        reaction = {}
                        for z in y:
                            if z.tag.split('}')[-1] == 'text':
                                reaction['Equation'] = z.text
                            elif 'type' in z.attrib and z.attrib['type'] == 'EC':
                                reaction['EC'] = z.attrib['id']
                        uniprot_data['Reaction'].append(reaction)

        # Process dbReference terms
        if 'dbReference' == x.tag.split('}')[-1]:

            # Process Go terms
            if 'type' in x.attrib and x.attrib['type'] == 'GO':
                uniprot_data['GO'][x.attrib['id']] = {}
                for y in x:
                    uniprot_data['GO'][x.attrib['id']][y.attrib['type']] = y.attrib['value']

            # Get PDB structures
            if 'type' in x.attrib and x.attrib['type'] == 'PDB':
                uniprot_data['Structures'].append(x.attrib['id'])

        if os.path.exists(xml_file):
            os.remove(xml_file)

    return uniprot_data

def retrieveAlphaFoldStructure(uniprot_code, output_directory, return_failed=False, version=4, index=None):
    """
    Download the AlphaFold structure for a specific Uniprot entry from the
    https://alphafold.ebi.ac.uk/ database.

    Parameters
    ==========
    uniprot_code : str
        Target uniprot code.
    output_directory : str
        Path to the output directory to store the file.
    index : str
        Index with ncbi code or alternative code, if any, to put in the output file name.
    Returns
    =======
    pdb_path : str
        Path to the alphafold structure.
    """

    # Create output directory
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # Remove extra info from up code
    if '_' in uniprot_code:
        uniprot_code_search = uniprot_code.split('_')[0]
    else:
        uniprot_code_search = uniprot_code

    # Define url address to retrieve alphafold structure
    url = 'https://alphafold.ebi.ac.uk/files/AF-'+uniprot_code_search+'-F1-model_v'+str(version)+'.pdb'

    # if index is passed, add the index to the filename
    if index:
        output_file = output_directory+'/AF-'+index+'-'+uniprot_code_search+'.pdb'
    else:
        output_file = output_directory+'/AF-'+uniprot_code_search+'.pdb'

    failed = False
    if not os.path.exists(output_file):
        try:
            print('Downloading Alphafold structure: AF-'+uniprot_code_search+'-F1-model_v'+str(version)+'.pdb')
            urllib.request.urlretrieve(url, output_file)
        except:
            print('Problem downloading AlphaFold structure for uniprot code %s' % uniprot_code)
            print('Please check that the following URL exists %s' % url)
            failed = True
            if return_failed:
                return failed

    else: # If file already exists
        print('Structure exists: '+output_file)

    if os.path.exists(output_file):
        if return_failed:
            return failed
        else:
            return output_file

class uniprot:

    def __init__(self, codes, data_folder='uniprot_data', database='uniprotkb', af_version=4):
        """
        Creates a class to download and parse uniprot data XML  and AlphaFold files.
        """

        # Check inputs
        if not isinstance(codes, (list, tuple, set)):
            raise ValueError('"codes" should be a list of Uniprot ID strings!')

        available_dbs = ['uniprotkb', 'uniref']
        if database not in available_dbs:
            raise ValueError('Wrong database selected. Available databases %s ' % available_dbs)

        # Set class attributes
        self.codes = codes
        self.data_folder = data_folder
        self.database = database
        self.xml_urls = {}
        self.xml_files = {}
        self.af_urls = {}
        self.af_files = {}
        self.fasta_urls = {}
        self.fasta_files = {}
        self.uniprot_data = {}
        self.sequences = {}

        # Create data folder
        if not os.path.exists(self.data_folder):
            os.mkdir(self.data_folder)

        # Create URL and XML paths
        for code in self.codes:
            self.xml_urls[code] = 'https://rest.uniprot.org/'+self.database+'/'+code+'.xml'
            self.xml_files[code] = self.data_folder+'/xml/'+code+'.xml'

        # Create URL and fasta files' paths
        for code in self.codes:
            self.fasta_urls[code] = 'https://rest.uniprot.org/'+self.database+'/'+code+'.fasta'
            self.fasta_files[code] = self.data_folder+'/fasta/'+code+'.fasta'

        # Create URL and AF models' paths
        for code in self.codes:
            self.af_urls[code] = 'https://alphafold.ebi.ac.uk/files/AF-'+code+'-F1-model_v'+str(af_version)+'.pdb'
            self.af_files[code] = self.data_folder+'/AF/'+code+'.pdb'

    def getXMLFiles(self):
        """
        Download uniprot XML files in parallel.
        """

        if not os.path.exists(self.data_folder+'/xml'):
            os.mkdir(self.data_folder+'/xml')

        jobs = [(self.xml_urls[code], self.xml_files[code]) for code in self.codes]
        pool = Pool(cpu_count())
        results = pool.map(self._download_file, jobs)

    def getFastaFiles(self):
        """
        Download uniprot fasta files in parallel.
        """

        if not os.path.exists(self.data_folder+'/fasta'):
            os.mkdir(self.data_folder+'/fasta')

        jobs = [(self.fasta_urls[code], self.fasta_files[code]) for code in self.codes]
        pool = Pool(cpu_count())
        results = pool.map(self._download_file, jobs)

    def getAlphaFoldFiles(self, version=None):
        """
        Download uniprot XML files in parallel.
        """

        if version != None:
            # Create URL and AF models' paths
            for code in self.codes:
                self.af_urls[code] = 'https://alphafold.ebi.ac.uk/files/AF-'+code+'-F1-model_v'+str(af_version)+'.pdb'
                self.af_files[code] = self.data_folder+'/AF/'+code+'.pdb'

        if not os.path.exists(self.data_folder+'/AF'):
            os.mkdir(self.data_folder+'/AF')

        jobs = [(self.af_urls[code], self.af_files[code]) for code in self.codes]
        pool = Pool(cpu_count())
        results = pool.map(self._download_file, jobs)

    def _download_file(self, url_and_file_path):
        """
        Download an url into the given file path

        Parameters
        ==========
        url : str
            The URL to use for downloading
        file_path : str
            Path where the file will be stored.
        """

        url = url_and_file_path[0]
        file_path = url_and_file_path[1]

        if os.path.exists(file_path):
            return True
        try:
            urllib.request.urlretrieve(url, file_path)
            return True
        except:
            print('Problem downloading file %s' % file_path)
            print('Please check that the following URL exists %s' % url)
            return False

    def parseXMLData(self):
        """
        Parse Uniprot XML files
        """

        for code in self.codes:
            if os.path.exists(self.xml_files[code]):
                self.uniprot_data[code] = self._parseXML(self.xml_files[code])
            else:
                self.uniprot_data[code] = None

        return self.uniprot_data

    def parseFastas(self):
        """
        Parse fasta files
        """

        for code in self.codes:
            if os.path.exists(self.fasta_files[code]):
                sequence = alignment.readFastaFile(self.fasta_files[code])
                for s in sequence:
                    self.sequences[code] = sequence[s]

        return self.sequences

    def _parseXML(self, xml_file):

        tree = ET.parse(xml_file)
        root = tree.getroot()

        # Create dictionary entries
        uniprot_data = {}
        uniprot_data['Protein'] = {}
        uniprot_data['Gene'] = {}
        uniprot_data['Organism'] = {}
        uniprot_data['GO'] = {}
        uniprot_data['Function'] = None
        uniprot_data['Reaction'] = []
        uniprot_data['Structures'] = []
        uniprot_data['features'] = {}
        tup = []

        # Iterate XML file to store UniProt data into a dictionary
        for x in root[0]:
            # Get sequence name
            if x.tag.endswith('name'):
                uniprot_data['Name'] = x.text

            # Get protein name information
            elif x.tag.endswith('protein'):
                for y in x:
                    if y.tag.endswith('recommendedName'):
                        for z in y:
                            if 'fullName' == z.tag.split('}')[-1]:
                                uniprot_data['Protein']['recommendedName'] = z.text
                    elif y.tag.endswith('alternativeName'):
                        for z in y:
                            if 'fullName' == z.tag.split('}')[-1]:
                                uniprot_data['Protein']['alternativeName'] = z.text

            # Get gene name information
            elif x.tag.endswith('gene'):
                for y in x:
                    if y.tag.endswith('name'):
                        uniprot_data['Gene'][y.attrib['type']] = y.text

            # Get organism information
            elif x.tag.endswith('organism'):
                uniprot_data['Organism']['Name'] = {}
                for y in x:
                    # Get organism names
                    if 'type' in y.attrib:
                        uniprot_data['Organism']['Name'][y.attrib['type']] = y.text
                    # Get taxonommic lineage
                    elif 'lineage' == y.tag.split('}')[-1]:
                        uniprot_data['Organism']['Lineage'] = []
                        for taxon in y:
                            uniprot_data['Organism']['Lineage'].append(taxon.text)

            #Get features information
            elif x.tag.endswith('feature'):
                t = x.attrib['type']
                if t not in uniprot_data['features']:
                    uniprot_data['features'][t] = []
                for y in x:
                    if y.tag.endswith('location'):
                        for p in y:
                            if p.tag.endswith('position'):
                                pos = p.attrib['position']
                                uniprot_data['features'][t].append(pos)
                            elif p.tag.startswith("begin") or p.tag.startswith("end"):
                                pos = p.attrib['position']
                                if len(tup) == 1:
                                    tup.append(pos)
                                    uniprot_data['features'][t].append(tup)
                                    tup = []
                                else:
                                    tup.append(pos)

            # Process comments
            elif 'comment' == x.tag.split('}')[-1]:

                # Get function
                if 'type' in x.attrib and x.attrib['type'] == 'function':
                    for y in x:
                        if 'text' == y.tag.split('}')[-1]:
                            uniprot_data['Function'] = y.text

                # Get catalytic activity
                if 'type' in x.attrib and x.attrib['type'] == 'catalytic activity':
                    for y in x:
                        if y.tag.split('}')[-1] == 'reaction':
                            reaction = {}
                            for z in y:
                                if z.tag.split('}')[-1] == 'text':
                                    reaction['Equation'] = z.text
                                elif 'type' in z.attrib and z.attrib['type'] == 'EC':
                                    reaction['EC'] = z.attrib['id']
                            uniprot_data['Reaction'].append(reaction)

            # Process dbReference terms
            if 'dbReference' == x.tag.split('}')[-1]:

                # Process Go terms
                if 'type' in x.attrib and x.attrib['type'] == 'GO':
                    uniprot_data['GO'][x.attrib['id']] = {}
                    for y in x:
                        uniprot_data['GO'][x.attrib['id']][y.attrib['type']] = y.attrib['value']

                # Get PDB structures
                if 'type' in x.attrib and x.attrib['type'] == 'PDB':
                    uniprot_data['Structures'].append(x.attrib['id'])

        return uniprot_data
