from HorusAPI import PluginPage, PluginEndpoint
from flask import request
import pele_analysis
import os

peleResultsPage = PluginPage(
    id="peleresults",
    name="Pele results",
    description="Pele results",
    html="peleResults.html"
)

def getPeleResults():
    try:
        data = request.json
    except Exception as e:
        print(e)

    if "peleSimulationFolder" in data:
        # TODO Create plots
        
        peleSimFolder = data['peleSimulationFolder']
        peleOutputFolder = data["peleOutputFolder"]

        if peleOutputFolder == '':
            peleOutputFolder = peleSimFolder + '../pele_output/'

        # add optinal value pele_data as output
        pele = pele_analysis.peleAnalysis(peleSimFolder, verbose=True, separator='-', trajectories=False, data_folder_name=peleOutputFolder, read_equilibration=True)

        print(pele)

        return {"ok": True, "msg": pele}

    return {"ok": False, "msg": "No pele simulations folder provided"}


# Add the endpoint to the PluginPage
testEndpoint = PluginEndpoint(
    url="/testEndpoint",
    methods=["POST"],
    function=getPeleResults,
)

# Add the endpoint to the page
peleResultsPage.addEndpoint(testEndpoint)