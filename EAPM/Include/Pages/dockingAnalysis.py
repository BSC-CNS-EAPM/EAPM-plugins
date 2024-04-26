from HorusAPI import PluginPage, PluginEndpoint
from flask import jsonify, request # type: ignore

docking_page = PluginPage(
    id="docking_analysis",
    name="Docking analysis",
    description="Analyze Glide results",
    html="docking_analysis.html",
)

def test():
    import math
    import json
    import pandas as pd

    import prepare_proteins
    import bsc_calculations

    from Bio import PDB

    data = request.json

    models_folder = data["modelsFolder"]
    docking_folder = data["dockingFolder"]

    models = prepare_proteins.proteinModels(models_folder)

    return jsonify({"dockingFolder": docking_folder, "modelsFolder": models_folder})

# Add the endpoint to the PluginPage
test_endpoint = PluginEndpoint(
    url="/api/test",
    methods=["POST"],
    function=test,
)

# Add the endpoint to the page
docking_page.addEndpoint(test_endpoint)