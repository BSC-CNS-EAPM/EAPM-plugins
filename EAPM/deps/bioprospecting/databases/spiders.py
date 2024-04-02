import scrapy
import json
import time

class PDBdatabase(scrapy.Spider):
    """
    Scrapy spider to retrieve information for a list of PDB ids. The PDB
    entry page is scrapped to extract information that is stored into a dictionary.
    The spider writes this dictionary into a json file.
    Attributes
    ----------
    pdb_ids : list
        List of uniprot ids to retrieve information.
    output_file : str
        Path for the output dictionary storing the retrieved information.
    """
    allowed_domains = ['www.rcsb.org/']

    def __init__(self, pdb_ids=None, output_file=None, **kwargs):
        self.pdb_ids = pdb_ids
        self.pdb_data = {}
        self.output_file = open(output_file, 'w')
        if self.pdb_ids == None:
            raise ValueError('You must give a list with the PDB IDs to retrieve\
                  information.')

    def start_requests(self):
        for pdbid in self.pdb_ids:
            yield scrapy.Request('https://www.rcsb.org/structure/'+pdbid, self.parse, meta={'pdbid': pdbid})

    def parse(self, response):

        # Get input uniprot id
        current = response.meta['pdbid']
        self.pdb_data[current] = {}

        # Save scraped url
        self.pdb_data[current]['url'] = 'https://www.rcsb.org/structure/'+current

        ### Scrape PDB data here ###

        ## Basic information
        self.parseBasicInformation(response, current)

        ## Macromolecules information
        self.parseMacromolecules(response, current)

        ## Small Molecules information
        self.parseSmallMolecules(response, current)

        ## Literature
        self.parseLiterature(response, current)

    def parseBasicInformation(self, response, current):
        ## Basic entry information ##
        structureTitle = response.css('#structureTitle::text').extract_first()
        self.pdb_data[current]['Title'] = structureTitle

        # Experimental data
        # for x in css.response('#exp_header_0_snapshot li').getall():
            # print(x)

        resolution = response.css('#exp_header_0_diffraction_resolution::text').extract_first()
        if resolution != None:
            resolution = float(resolution.replace('Å',''))
        self.pdb_data[current]['Resolution'] = resolution

    def parseMacromolecules(self, response, current):
        ## Macromolecules entry information ##
        self.pdb_data[current]['Macromolecules'] = {}
        molecule_names = []
        chains = []
        organisms = []
        lengths = []
        uniprot_ids = []

        # Scrape macromolecules table
        for i,x in enumerate(response.css('#MacromoleculeTable tr:nth-child(3) td')):
            index = (i+1)%6
            if index== 1:
                molecule_name = x.css('td::text').extract_first()
                molecule_names.append(molecule_name)
            elif index == 2:
                chain_ids = x.css('a::text').extract()
                chains.append(chain_ids)
            elif index == 3:
                length = x.css('td::text').extract_first()
                lengths.append(length)
            elif index == 4:
                organism = x.css('td a::text').extract_first()
                organisms.append(organism)

        for i,x in enumerate(response.css('table')):
            id_data = x.css('::attr(id)').extract_first()
            if id_data != None and 'table_macromolecule-protein-entityId' in id_data:
                upid = x.css('.text-left .querySearchLink::text').extract()
                if upid == []:
                    uniprot_ids.append(None)
                else:
                    uniprot_ids.append(upid[0])

        # Store data in dictionary
        for entity in zip(chains, molecule_names, lengths, organisms, uniprot_ids):
            for chain in entity[0]:
                self.pdb_data[current]['Macromolecules'][chain] = {}
                self.pdb_data[current]['Macromolecules'][chain]['Molecule'] = entity[1]
                self.pdb_data[current]['Macromolecules'][chain]['Length'] = entity[2]
                self.pdb_data[current]['Macromolecules'][chain]['Organism'] = entity[3]
                self.pdb_data[current]['Macromolecules'][chain]['Uniprot'] = entity[4]

    def parseSmallMolecules(self, response, current):
        ## Small Molecules entry information ##
        self.pdb_data[current]['Small Molecules'] = {}
        ligand_ids = []
        chains = []
        names = []
        for i,x in enumerate(response.css('#LigandsMainTable td')):
            if (i+1)%5 == 1:
                ligand_ids.append(x.css('a::text').extract_first())
            elif (i+1)%5 == 2:
                chains.append(x.css('div.ellipsisToolTip::text').extract())
            elif (i+1)%5 == 3:
                names.append(x.css('strong::text').extract_first())

        for chain in chains:
            for c in chain:
                self.pdb_data[current]['Small Molecules'][c] = {}
                self.pdb_data[current]['Small Molecules'][c]['ID'] = []
                self.pdb_data[current]['Small Molecules'][c]['Name'] = []

        for entity in zip(chains, ligand_ids, names):
            for chain in entity[0]:
                self.pdb_data[current]['Small Molecules'][chain]['ID'].append(entity[1])
                self.pdb_data[current]['Small Molecules'][chain]['Name'].append(entity[2])

    def parseLiterature(self, response, current):
        try:
            self.pdb_data[current]['Literature'] = {}
            litpanel = response.css('#literaturepanel .panel-body')
            title = litpanel.css('#primarycitation h4::text').extract()[0]
            self.pdb_data[current]['Literature']['Title'] = title
            for a in litpanel.css('#primarycitation a'):
                href = a.css('::attr(href)').extract_first()
                if 'dx.doi.org' in href:
                    self.pdb_data[current]['Literature']['DOI'] = href
                    break
        except:
            self.pdb_data[current]['Literature'] = None

    def closed(self, spider):
        json.dump(self.pdb_data, self.output_file)
        self.output_file.close()
