from HorusAPI import PluginPage, PluginEndpoint
from flask import request
import os

peleResultsPage = PluginPage(
    id="peleresults",
    name="Pele results",
    description="Pele results",
    html="peleResults.html"
)





def getPeleResults():
    import pele_analysis
    from flask import jsonify
    try:
        data = request.json
    

        if "peleSimulationFolder" in data:
            # TODO Create plots
            
            peleSimFolder = data['peleSimulationFolder']

            if not peleSimFolder.endswith('/'):
                peleSimFolder += '/'

            # peleOutputFolder = data["peleOutputFolder"]

            # if peleOutputFolder == '':
            #     peleOutputFolder = peleSimFolder + '../pele_output/'

            peleOutputFolder = 'pele_output'

            containerFoler = '/'.join(peleSimFolder.split('/')[:-2])

            os.chdir(containerFoler)

            peleFolder = peleSimFolder.split("/")[len(peleSimFolder.split("/")) - 2]
            # pele = pele_analysis.peleAnalysis(peleFolder, verbose=True, separator='-', trajectories=False, data_folder_name=peleOutputFolder, read_equilibration=True)
            pele = pele_analysis.peleAnalysis('pele/', verbose=True, separator='-', trajectories=False, data_folder_name='pele_data/', read_equilibration=True)
            print('saliendo analysis')

            # Classify distances into common metrics --> the catalytic labels
            catalytic_names = ['SER-L', 'SER-HIS', 'HIS-ASP']

            catalytic_labels = {}
            for cn in catalytic_names:
                catalytic_labels[cn] = {}
                for protein in pele.proteins:
                    if protein not in catalytic_labels[cn]:
                        catalytic_labels[cn][protein] = {}
                    for ligand in pele.ligands:
                        if ligand not in catalytic_labels[cn][protein]:
                            catalytic_labels[cn][protein][ligand] = []

            for protein in pele.proteins:
                for ligand in pele.ligands:
                    distances = pele.getDistances(protein, ligand)
                    if distances == None:
                        continue
                    for d in distances:
                        at1 = d.replace('distance_','').split('_')[0]
                        at2 = d.replace('distance_','').split('_')[1]
                        if at1.endswith('OG') and at2.startswith('L'):
                            catalytic_labels['SER-L'][protein][ligand].append(d)

                        elif at1.endswith('OG') and at2.endswith('NE2'):
                            catalytic_labels['SER-HIS'][protein][ligand].append(d)

                        elif at1.endswith('ND1') and (at2.endswith('OD1') or at2.endswith('OD2')):
                            catalytic_labels['HIS-ASP'][protein][ligand].append(d)

            pele.combineDistancesIntoMetrics(catalytic_labels, overwrite=True)

            return jsonify({"ok": True, "distances": distances, "proteins": pele.proteins, "ligands": pele.ligands})

        return jsonify({"ok": False, "msg": "No pele simulations folder provided"})
    except Exception as e:
            return jsonify({"ok": False, "msg": str(e)})

# Add the endpoint to the PluginPage
testEndpoint = PluginEndpoint(
    url="/testEndpoint",
    methods=["POST"],
    function=getPeleResults,
)

# Add the endpoint to the page
peleResultsPage.addEndpoint(testEndpoint)