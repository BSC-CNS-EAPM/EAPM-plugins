import os
import numpy as np
from . import alignment
from . import databases
from . import ssn

class protocol:

    def __init__(self, target_sequences):

        self.target_sequences = target_sequences
        self.directories = {}
        self.directories['base'] = '.bioprospecting'
        self.directories['psiblast_results'] = {}
        self.psiblast_results = {}
        self.databases = ['Uniprot']
        self.ssn_matrix_base = {}
        self.sequences = None
        self.cluster_score = None

        for database in self.databases:
            self.directories['psiblast_results'][database] = self.directories['base']+'/psiblast_results/'+database

        # Create bioprosepcting folders
        if not os.path.exists(self.directories['base']):
            os.mkdir(self.directories['base'])

        if not os.path.exists(self.directories['base']+'/psiblast_results'):
            os.mkdir(self.directories['base']+'/psiblast_results')

        for database in self.databases:
            if not os.path.exists(self.directories['base']+'/psiblast_results/'+database):
                os.mkdir(self.directories['base']+'/psiblast_results/'+database)

        if not os.path.exists(self.directories['base']+'/psiblast_results/'):
            os.mkdir(self.directories['base']+'/psiblast_results/')

        # PSI blast databases
        for database  in self.databases:

            self.psiblast_results[database] = {}

            for code in self.target_sequences:
                if database == 'Uniprot':
                    pbr_dir = self.directories['psiblast_results'][database]+'/'+code+'.json'
                    if not os.path.exists(pbr_dir):
                        self.psiblast_results['Uniprot'][code] = databases.PSIBlastUniProtDatabase(target_sequences[code])
                        alignment.savePSIBlastAsJson(self.psiblast_results['Uniprot'][code],
                                                     pbr_dir)
                    else:
                        self.psiblast_results['Uniprot'][code] = alignment.readPSIBlastFromJson(pbr_dir)

    def calculateSSN(self, prot_databases=None, overwrite=False):

        if prot_databases == None:
            prot_databases = self.databases

        self.sequences = {}

        all_codes = []
        for database in prot_databases:
            for code in self.psiblast_results[database]:
                for i in self.psiblast_results[database][code]:
                    for pb_code in self.psiblast_results[database][code][i]:
                        #all_codes.append(pb_code)
                        all_codes.append(pb_code.split('|')[1])


        self.sequences = databases.getUniprotSequences(all_codes)

        for code,sequence in self.target_sequences.items():
            self.sequences[code] = sequence

        self.ssn_matrix_base = ssn.sequenceSimilarityNetwork(self.sequences,
                                                        similarity_matrix_file=self.directories['base']+'/ssn_matrix.npy',
                                                        overwrite=overwrite, target_sequences = list(self.target_sequences.keys()))

    def cluster_threshold_analysis(self, threshold_min, threshold_max, attribute ,step=0.01, display=True):
        threshold_list = np.around(np.arange(threshold_min,threshold_max+step,step),decimals = 2)
        self.ssn_matrix_base.getUniProtAttributes()
        self.ssn_matrix_base.colorNodeFillByAttribute(attribute,overwrite=True)
        self.cluster_score = self.ssn_matrix_base.clusterAnalysis(threshold_list,attribute)

        if display == True:
            self.ssn_matrix_base.createNetwork(threshold_list)
            self.ssn_matrix_base.cleanDisplay(4)
            self.ssn_matrix_base.drawInteractiveNetwork(min(self.cluster_score,key=self.cluster_score.get),threshold_min,threshold_max,step,score = self.cluster_score)
