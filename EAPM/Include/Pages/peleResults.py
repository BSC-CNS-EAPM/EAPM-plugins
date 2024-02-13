from HorusAPI import PluginPage, PluginEndpoint
from flask import request
import os

peleResultsPage = PluginPage(
    id="peleresults",
    name="Pele results",
    description="Pele results",
    html="peleResults.html"
)

def returnDataframe(self, protein, ligand, vertical_line=None, color_column=None, size=1.0, labels_size=10.0, plot_label=None,
                                        xlim=None, ylim=None, metrics=None, labels=None, title=None, title_size=14.0, return_axis=False, dpi=300, show_legend=False,
                                        axis=None, xlabel=None, ylabel=None, vertical_line_color='k', vertical_line_width=0.5, marker_size=0.8, clim=None, show=False,
                                        clabel=None, legend_font_size=6, no_xticks=False, no_yticks=False, no_cbar=False, no_xlabel=False, no_ylabel=False,
                                        relative_color_values=False, dataframe=None, **kwargs):
        """
        Creates a scatter plot for the selected protein and ligand using the x and y
        columns. Data series can be filtered by specific metrics.

        Parameters
        ==========
        protein : str
            The target protein.
        ligand : str
            The target ligand.
        x : str
            The column name of the data to plot in the x-axis.
        y : str
            The column name of the data to plot in the y-axis.
        vertical_line : float
            Position to plot a vertical line.
        color_column : str
            The column name to use for coloring the plot. Also a color cna be given
            to use uniformly for the points.
        xlim : tuple
            The limits for the x-range.
        ylim : tuple
            The limits for the y-range.
        clim : tuple
            The limits for the color range.
        metrics : dict
            A set of metrics for filtering the data points.
        labels : dict
            Analog to metrics, use the label column values to filter the data.
        title : str
            The plot's title.
        return_axis : bool
            Whether to return the axis of this plot.
        axis : matplotlib.pyplot.axis
            The axis to use for plotting the data.
        """
        try:
            if not isinstance(dataframe, type(None)):
                protein_series = dataframe[dataframe.index.get_level_values('Protein') == protein]
            else:
                protein_series = self.data[self.data.index.get_level_values('Protein') == protein]

            if protein_series.empty:
                raise ValueError('Protein name %s not found in data!' % protein)
            ligand_series = protein_series[protein_series.index.get_level_values('Ligand') == ligand]

            if ligand_series.empty:
                raise ValueError("Ligand name %s not found in protein's %s data!" % (ligand, protein))

            # Add distance data to ligand_series
            if len(ligand_series) != 0:
                if protein in self.distances:
                    if ligand in self.distances[protein]:
                        if not isinstance(self.distances[protein][ligand], type(None)):
                            for distance in self.distances[protein][ligand]:
                                #if distance.startswith('distance_'):

                                if not isinstance(dataframe, type(None)):
                                    
                                    indexes = dataframe.reset_index().set_index(['Protein', 'Ligand', 'Epoch', 'Trajectory', 'Accepted Pele Steps', 'Step']).index
                                    ligand_series[distance] = self.distances[protein][ligand][self.distances[protein][ligand].index.isin(indexes)][distance].tolist()
                                else:
                                    ligand_series[distance] = self.distances[protein][ligand][distance].tolist()

            # Filter points by metric
            if not isinstance(metrics, type(None)):
                for metric in metrics:
                    mask = ligand_series[metric] <= metrics[metric]
                    ligand_series = ligand_series[mask]

            if not isinstance(labels, type(None)):
                for label in labels:
                    if labels[label] != None:
                        mask = ligand_series[label] == labels[label]
                        ligand_series = ligand_series[mask]
            return ligand_series
        except Exception as e:
            print(e)
            return {"error": e}

def getPeleResults():
    import pele_analysis
    from flask import jsonify

    try:
        data = request.json

        if "peleSimulationFolder" in data:
            if data['desiredProtein'] == None:
                peleSimFolder = data['peleSimulationFolder']

                if not peleSimFolder.endswith('/'):
                    peleSimFolder += '/'

                peleOutputFolder = data["peleOutputFolder"]

                if peleOutputFolder == '':
                    peleOutputFolder = peleSimFolder + '../pele_data/'

                containerFoler = '/'.join(peleSimFolder.split('/')[:-2])

                os.chdir(containerFoler)

                peleFolder = peleSimFolder.split("/")[len(peleSimFolder.split("/")) - 2]
                pele = pele_analysis.peleAnalysis(peleFolder, verbose=True, separator='-', trajectories=False, data_folder_name=peleOutputFolder, read_equilibration=True)
                
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


                distances = {}
                for protein in pele.proteins:
                    distances[protein] = {}
                    for ligand in pele.ligands:
                        gettedDistances = pele.getDistances(protein, ligand)
                        if gettedDistances != []:
                            distances[protein][ligand] = gettedDistances

                firstProteinCouple = None
                firstLigandCouple = None

                for protein in pele.proteins:
                    for ligand in pele.ligands:
                        if pele.getDistances(protein, ligand):
                            firstProteinCouple = protein
                            firstLigandCouple = ligand
                            break
                    if firstProteinCouple:
                        break

                pele.combineDistancesIntoMetrics(catalytic_labels, overwrite=True)

                dict = returnDataframe(pele, firstProteinCouple, firstLigandCouple)

                df_flat_index = dict.reset_index()

                df_flat_index.fillna(0, inplace=True)

                data_dict = df_flat_index.to_dict(orient='split')

                return jsonify({"ok": True, "distances": distances, "proteins": pele.proteins, "ligands": pele.ligands, "dict": data_dict})
            else:
                peleSimFolder = data['peleSimulationFolder']

                if not peleSimFolder.endswith('/'):
                    peleSimFolder += '/'

                peleOutputFolder = data["peleOutputFolder"]

                if peleOutputFolder == '':
                    peleOutputFolder = peleSimFolder + '../pele_data/'

                containerFoler = '/'.join(peleSimFolder.split('/')[:-2])

                os.chdir(containerFoler)

                peleFolder = peleSimFolder.split("/")[len(peleSimFolder.split("/")) - 2]
                pele = pele_analysis.peleAnalysis(peleFolder, verbose=True, separator='-', trajectories=False, data_folder_name=peleOutputFolder, read_equilibration=True)

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

                distances = {}
                for protein in pele.proteins:
                    distances[protein] = {}
                    for ligand in pele.ligands:
                        gettedDistances = pele.getDistances(protein, ligand)
                        if gettedDistances != []:
                            distances[protein][ligand] = gettedDistances

                # TODO IndexOf and get the desired
                            
                desiredProtein = data['desiredProtein']
                desiredLigand = data['desiredLigand']

                pele.combineDistancesIntoMetrics(catalytic_labels, overwrite=True)

                dict = returnDataframe(pele, pele.proteins[pele.proteins.index(desiredProtein)], pele.ligands[pele.ligands.index(desiredLigand)])

                df_flat_index = dict.reset_index()

                data_dict = df_flat_index.to_dict(orient='split')

                return jsonify({"ok": True, "dict": data_dict, "distances": distances})

        return jsonify({"ok": False, "msg": "No pele simulations folder provided"})
    except Exception as e:
            print(f'Error: {e}')
            return jsonify({"ok": False, "msg": str(e)})

def getPlotlyCLick():
    import pele_analysis
    from flask import jsonify

    data = request.json

    if "distance" and "bindingEnergy" in data:

        distance = data['distance']
        bindingEnergy = data['bindingEnergy']

        print(f'{distance} - {bindingEnergy}')
        return jsonify({"ok": True, "msg": f'{distance} - {bindingEnergy}'})

    try:
        pass
    except Exception as e:
        print(f'Error: {e}')
        return jsonify({"ok": False, "msg": str(e)})

# Add the endpoint to the PluginPage
peleResults = PluginEndpoint(
    url="/peleResults",
    methods=["POST"],
    function=getPeleResults,
)

plotlyCLick = PluginEndpoint(
    url="/plotlyClick",
    methods=["POST"],
    function=getPlotlyCLick,
)

# Add the endpoint to the page
peleResultsPage.addEndpoint(peleResults)
peleResultsPage.addEndpoint(plotlyCLick)
