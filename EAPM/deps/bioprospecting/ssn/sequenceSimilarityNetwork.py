import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import ipywidgets as widgets
from ipywidgets import interact
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage, fcluster
from collections import Counter
from math import log as ln

from .. import alignment
from .. import databases

class sequenceSimilarityNetwork:
    """
    Class to perform sequence similarity networks.

    Attributes
    ----------
    sequences : dict
        Input dictionary mapping each entry code to its sequence.
    target_sequences : array
        List of the codes of the target sequences.
    codes : list
        Reference list containing the protein codes in the order they are in the
        similarity matrix.
    similarity_matrix : np.ndarray
        Array containing the similarity values between each pair of sequences.
    node_attributes
        Dictionary holding attributes add to each node.
    networks
        Dictionary holding networks generated for each given threshold.
    network_display
        Dictionary of drawn network objects.
    color_dictionary
        Dictionary of attribute for each color

    Methods
    -------
    plotSimilarityMatrix()
        Generates a plot of the similarity matrix.
    createNetwork()
        Creates a network for a specific threshold and adds it to the networks dictionary
    createSifNetworkFile()
        Creates a network sif format file for a specific threshold.
    addNodeAttribute()
        Add attributes to the nodes in the network.
    createNodeAttributesFile()
        Generates a file containing all the attributes for the nodes in the network.
    """

    def __init__(self, sequences, similarity_matrix=None, similarity_matrix_file='pid.matrix.npy',
                 write_pid_file=True, overwrite=False, verbose=False,target_sequences = []):
        """
        Initializes the sequenceSimilarityNetwork class by giving a dictionary
        containing the protein codes as keys and sequences as values. The methods
        calculates and stores a PID matrix in a numpy format file. If the pid_file
        path is found then the matrix is read from this file, otherwise the PID
        matrix is calculated and the content is stored in the pid_file path. Alternatively
        a precalcualted similarity matrix can be given in numpy format.

        Parameters
        ----------
        sequences : dict
            Dictionary mapping each entry protein code to its sequence.
        similarity_matrix : numpy.ndarray
            The similarity matrix for the input sequences as a numpy array.
        similarity_matrix_file : str
            Path to the similarity matrix file in numpy format.
        write_pid_file : bool
            Write the PID matrix to a file?
        overwrite : bool
            Forces the recalcualation of the PID matrix and overwrites the pid_file
            file.
        """

        # Create class variables
        self.sequences = sequences
        if not similarity_matrix_file.endswith('.npy'):
            similarity_matrix_file = similarity_matrix_file+'.npy'
        self.target_sequences = target_sequences
        self.similarity_matrix_file = similarity_matrix_file
        self.node_attributes = {}
        self.networks = {}
        self.networks_display = {}
        self.color_dictionary = {}
        # Store sequences and order of codes
        self.codes = [code for code in self.sequences]

        # Check input parameters
        if not isinstance(sequences, dict):
            raise ValueError('Input sequences should be a dictionary.')

        if similarity_matrix != None:
            self.similarity_matrix = similarity_matrix
            overwrite = True # Write the similarity matrix into a file
        else:
            # Check if the similarity matrix file exists
            if os.path.exists(self.similarity_matrix_file) and not overwrite:
                print('Similarity matrix file found! Reading the similarity matrix from '+
                self.similarity_matrix_file)
                self.similarity_matrix = np.load(self.similarity_matrix_file)
            else:
                # Calcualte the similarity matrix
                self.similarity_matrix = np.zeros((len(self.codes), len(self.codes)))
                # Calculate similarity matrix using blast
                sequences = [self.sequences[code] for code in self.codes]
                for i in range(len(self.codes)):
                    if verbose:
                        print('Calculating PID values for : '+self.codes[i])
                    self.similarity_matrix[i] = alignment.blast.calculatePIDs(sequences[i],
                                                sequences)
                overwrite = True # Write the similarity matrix into a file

        # Check problems with the similarity matrix
        if self.similarity_matrix.shape[0] != self.similarity_matrix.shape[1]:
            raise ValueError('The matrix in the file is not a square matrix')
        if self.similarity_matrix.shape[0] != len(self.codes):
            raise ValueError('The first matrix dimmension do not mathch the number of sequences')
        if self.similarity_matrix.shape[1] != len(self.codes):
            raise ValueError('The seconde matrix dimmension do not mathch the number of sequences')

        if overwrite:
            np.save(self.similarity_matrix_file, self.similarity_matrix)

    def createNetwork(self, threshold_list):
        """
        Creates networks for the given thresholds and adds them to the networks dictionary

        Parameters
        ----------
        threshold_list : list
            List of similarity values to create the edges of each network.
        """

        for threshold in threshold_list:

            if threshold < 0 or threshold > 1:
                raise ValueError('Wrong Threshold! The threshold must be a number between zero and one.')

            G = nx.Graph()
            G.add_nodes_from(self.codes)
            nx.set_node_attributes(G,self.node_attributes)

            added = []
            # Write edges above threshold
            for i in range(self.similarity_matrix.shape[0]):
                for j in range(self.similarity_matrix.shape[1]):
                    if j > i:
                        if self.similarity_matrix[i][j] >= threshold:
                            G.add_edge(self.codes[i],self.codes[j])
                            if self.codes[i] not in added:
                                added.append(self.codes[i])
                            if self.codes[j] not in added:
                                added.append(self.codes[j])

            self.networks[threshold] = G
            self.networks_display[threshold] = nx.spring_layout(G,k=None,iterations=50)

    def createSifNetworkFile(self, threshold, output_file='network.sif', overwrite=None):
        """
        Generates a network file based on the given threshold.

        Parameters
        ----------
        threshold : float
            Similarity value to create the network edges.
        output_file : str
            Name of the network output file.
        overwrite : bool
            Whether to overwrite an existing network file.
        """

        if threshold < 0 or threshold > 1:
            raise ValueError('Wrong Threshold! The threshold must be a number between zero and one.')

        if not output_file.endswith('.sif'):
            output_file = output_file+'.sif'

        if output_file == 'network.sif':
            output_file = output_file.replace('.sif', str(threshold)+'.sif')

        if os.path.exists(output_file) and not overwrite:
            raise ValueError('Network file exists. Give overwrite=True to replace it.')


        added = []
        with open(output_file, 'w') as nf:
            # Write edges above threshold
            for i in range(self.similarity_matrix.shape[0]):
                for j in range(self.similarity_matrix.shape[1]):
                    if j > i:
                        if self.similarity_matrix[i][j] >= threshold:
                            nf.write('%s 1.0 %s\n' % (self.codes[i], self.codes[j]))
                            if self.codes[i] not in added:
                                added.append(self.codes[i])
                            if self.codes[j] not in added:
                                added.append(self.codes[j])

            # Write remaining codes
            for code in self.codes:
                if code not in added:
                    nf.write('%s\n' % code)


    def addNodeAttribute(self, attribute_name, attribute_values, overwrite=False):
        """
        Add an attribute category to the proteins in the network. The input is a
        dictionary mapping the attribute values to the protein codes. All missing
        codes will be given a 'none' value. If the attribute name or any attribute
        value have a space " " character it will be replaced by an underscore "_"
        character.

        Parameters
        ----------
        attribute_name : str
            Name of the node attribute to be added.
        attribute_values : dict
            Dictionary mapping the protein codes to the attribute values.
        overwrite : bool
            Whether to overwrite the previous attribute values with the same attributes
            name.
        """

        # Remove empty spaces from the attribute name

        if ' ' in attribute_name:
            print('Space " " found in attribute %s' % attribute_name)
            attribute_name = attribute_name.replace(' ', '_')
            print('The new attribute name will be %s' % attribute_name)


        # Check for nodes not present in the network.
        for code in attribute_values:
            if code not in self.codes:
                print('Code %s is not present in the network!' % code)

        # Iterate for codes
        for code in self.codes:
            if code not in self.node_attributes:
                self.node_attributes[code] = {}

            # check if attribute already exists
            if attribute_name in self.node_attributes[code] and not overwrite:
                raise ValueError('Attribute name already found. Give overwrite=True to rewrite the values.')

            if code not in attribute_values:
                self.node_attributes[code][attribute_name] = 'none'
            else:
                if isinstance(attribute_values, str):
                    if ' ' in attribute_values[code]:
                        print('Space " " found in attribute value %s' % attribute_values[code])
                        attribute_values[code] = attribute_values[code].replace(' ', '_')
                        print('The new attribute value will be %s' % attribute_values[code])
                    self.node_attributes[code][attribute_name] = attribute_values[code]
                else:
                    self.node_attributes[code][attribute_name] = str(attribute_values[code])

    def addSequenceLengthAttribute(self):
        """
        Add the length of each sequence as a node attribute in the network.
        """

        # Iterate for codes
        for code in self.codes:
            if code not in self.node_attributes:
                self.node_attributes[code] = {}
            self.node_attributes[code]['sequence_length'] = str(len(self.sequences[code]))
            #self.node_attributes[code]['sequence_length'] = str(len(self.sequences[code]))[0] + "00"

    def getUniProtAttributes(self):

        attribute_dict = {}

        for code in self.codes:
            if code not in self.target_sequences:
                attributes = (databases.getUniprotData(code))
                for att_name,att_value in attributes.items():
                    if att_name not in attribute_dict:
                        attribute_dict[att_name] = {}
                    attribute_dict[att_name][code] = att_value

        for key,value in attribute_dict.items():
            self.addNodeAttribute(key,value,overwrite=True)

    def createNodeAttributesFile(self, output_file, overwrite=None, separator=', '):
        """
        Create a node property file using the properties given to the sequenceSimilarityNetwork.

        Parameters
        ----------
        output_file : str
            Name of the node property output file.
        overwrite : bool
            Whether to overwrite an existing node attributes file.
        """

        if os.path.exists(output_file) and not overwrite:
            raise ValueError('Attribute file exists. Give overwrite=True to replace it.')

        with open(output_file, 'w') as af:

            # Get attributes list
            attributes = []
            for code in self.node_attributes:
                for attribute in self.node_attributes[code]:
                    attributes.append(attribute)
                break

            # Write header line
            af.write('code'+separator+separator.join(attributes)+'\n')

            for code in self.codes:
                attribute_line = ''
                for i,attribute in enumerate(attributes):
                    attribute_line += self.node_attributes[code][attribute].replace(separator.replace(' ',''), '.')
                    if i != len(attributes)-1:
                        attribute_line += separator
                af.write(code+separator+attribute_line+'\n')

    def colorNodeFillByAttribute(self, attribute_name, overwrite=False):
        """
        Create a node attribute to color nodes based on the value of an attribute.
        The colors are limited to 10 different attribute names, they are colored from
        the most numerous to the less numerous. Further than 10 attribute names will
        not be assigned a color.

        Parameters
        ----------
        attribute_name : str
            Name of the attribute to use to assign colors
        """

        color_palette = [matplotlib.colors.rgb2hex(matplotlib.pyplot.cm.tab10(i)) for i in range(10)]

        if self.node_attributes == {}:
            raise ValueError('There are no attributes assigned. Use addNodeAttribute to add attributes.')

        names = []
        for code in self.codes:
            if attribute_name in self.node_attributes[code]:
                names.append(self.node_attributes[code][attribute_name])
            else:
                raise ValueError('attribute_name %s not found in node attributes' % attribute_name)

        #names = [self.node_attributes[attribute_name][n] for n in self.node_attributes[attribute_name]]
        names_set = set(names)
        names_count = sorted([ (name,names.count(name)) for name in names_set], key=lambda x:x[1], reverse=True)
        colors = {}

        for i,name in enumerate(names_count):

            for code in self.codes:
                if self.node_attributes[code][attribute_name] == name[0]:
                    if i < 9:
                        colors[code] = color_palette[i]
                        self.color_dictionary[color_palette[i]] = self.node_attributes[code][attribute_name][0:50]
                    elif code in self.target_sequences:
                        colors[code] = matplotlib.colors.rgb2hex(matplotlib.pyplot.cm.tab10(10))
                        self.color_dictionary[matplotlib.colors.rgb2hex(matplotlib.pyplot.cm.tab10(10))] = "Target"
                    else:
                        colors[code] = '#000000'
                        self.color_dictionary['#000000'] = "Others"

        #print(colors)

        self.addNodeAttribute('fill_color', colors, overwrite=overwrite)


    def cleanDisplay(self,min_comp_size=3):
        """
        Clean the display moving all unconnected nodes to the side
        """

        for threshold in self.networks.keys():

            uncon = list([list(comp) for comp in nx.connected_components(self.networks[threshold]) if len(comp) < min_comp_size])
            uncon_nodes = np.array([code for code,coord in self.networks_display[threshold].items() if (any(code in comp for comp in uncon))])
            con_nodes = np.array([code for code,coord in self.networks_display[threshold].items() if not(any(code in comp for comp in uncon))])

            con_graph = self.networks[threshold].copy()
            for node in uncon_nodes:
                con_graph.remove_node(node)

            pos = nx.spring_layout(con_graph,k=None,iterations=50)

            for node in con_nodes:
                self.networks_display[threshold][node] = pos[node]

            coordinates = np.array([coord for code,coord in pos.items()])

            xmax, ymax = coordinates.max(axis=0)
            xmin, ymin = coordinates.min(axis=0)
            x_cord = xmax +0.05
            y_cord = ymax

            rearanged = []
            for comp in uncon:
                for i,node in enumerate(comp):
                    if i == 0:
                        ref_position = self.networks_display[threshold][node]
                        self.networks_display[threshold][node] = np.array([x_cord,y_cord])
                    else:
                        mov_x = x_cord - ref_position[0]
                        mov_y = y_cord - ref_position[1]
                        self.networks_display[threshold][node] += np.array([mov_x,mov_y])


                if y_cord < ymin:
                    x_cord += 0.1
                    y_cord = ymax
                else:
                    y_cord -= 0.1



    #def drawNetwork(self,network_thresholds = 'all'):
    def drawNetwork(self,threshold):
        """
        Draw the networks for each threshold passed. The layout of the network is based on a a force-directed algorithm

        Parameters
        ----------
        threshold_list : list
            List of similarity values of the created networks. If no argument is given all networks are drawn.
        """

        #if network_thresholds == 'all':
        #    network_thresholds = self.networks.keys()
        #else:
        #    if any(item not in self.networks.keys() for item in network_thresholds):
        #        raise ValueError('Network with threshold specified in list has not been created. You can use createNetwork function to do it.')


        #for threshold,network in self.networks.items():
        #    if threshold in network_thresholds:
        #        color_map = nx.get_node_attributes(network,"fill_color")
        #        pos = nx.spring_layout(network,k=None,iterations=50)
        #        fig = plt.figure(figsize=(30, 30),dpi=200)
        #        nx.draw_networkx(network,pos=pos,node_color=color_map.values(), with_labels=False,node_size = 50,width = 0.25)
        #        #nx.draw_spring(network,node_color=color_map.values(), with_labels=False)
        #        plt.title(str(threshold))
        #        plt.show()


        if threshold not in self.networks:
             raise ValueError('Network with threshold %s has not been created. You can use createNetwork function to do it.' % str(threshold))


        network = self.networks[threshold]

        color_map = nx.get_node_attributes(network,"fill_color")
        pos = nx.spring_layout(network,k=None,iterations=50)
        self.networks_display[threshold] = nx.spring_layout(network,k=None,iterations=50)

        fig = plt.figure(figsize=(30, 30),dpi=200)
        nx.draw_networkx(network,pos=pos,node_color=color_map.values(), with_labels=False,node_size = 50,width = 0.25)
        #nx.draw_spring(network,node_color=color_map.values(), with_labels=False)
        plt.title(str(threshold))
        plt.show()

    def drawInteractiveNetwork(self,threshold,min,max,step,score=None):

        if threshold not in self.networks:
            raise ValueError('Network with threshold %s has not been created. You can use createNetwork function to do it.' % str(threshold))

        def draw(threshold):
            color_map = nx.get_node_attributes(self.networks[threshold],"fill_color")
            #fig = plt.figure(figsize=(15, 15),dpi=100)
            #nx.draw_networkx(self.networks[threshold],pos=self.networks_display[threshold],node_color=color_map.values(), with_labels=False,node_size = 50,width = 0.25)
            #nx.draw_spring(network,node_color=color_map.values(), with_labels=False)
            fig, ax = plt.subplots(1, figsize=(15, 15),dpi=100)

            color_node = {}
            for key, value in color_map.items():
                if value in color_node:
                    color_node[value].append(key)
                else:
                    color_node[value]=[key]

            for color,nodes in color_node.items():
                nlist = nodes
                ncolor = color

                #for node in nodes:
                #    if node in self.target_sequences:
                #        nx.draw_networkx_nodes(self.networks[threshold],pos = self.networks_display[threshold],nodelist=[node],node_color=ncolor,ax=ax,label=self.color_dictionary[color],node_size = 100)
                #        nlist.remove(node)

                nx.draw_networkx_nodes(self.networks[threshold],pos = self.networks_display[threshold],nodelist=nlist,node_color=ncolor,ax=ax,label=self.color_dictionary[color],node_size = 50)


            nx.draw_networkx_edges(self.networks[threshold],pos = self.networks_display[threshold],ax=ax,width = 0.25)
            ax.legend(scatterpoints=1)

            if score == None:
                plt.title(str(threshold))
            else:
                plt.title("Threshold: " + str(threshold) + ", Score: " + str(score[threshold]))
            plt.show()

        interact(draw,threshold=widgets.FloatSlider(min=min, max=max, step=step, value=threshold))


    def clusterAnalysis(self,threshold_list,attribute):

        results = {}
        dist_matrix = (1-self.similarity_matrix)
        y = ssd.squareform(dist_matrix,checks=False)
        Z = linkage(y,method='ward', optimal_ordering=False)

        for threshold in threshold_list:
            clusters = fcluster(Z,threshold, criterion='distance')
            clusters = list(map(str, clusters))
            zip_iterator = zip(self.codes, (clusters))
            clus_dic = dict(zip_iterator)

            self.addNodeAttribute("clusters",clus_dic,overwrite=True)

            att_frequency_clus = {}
            counts = Counter(clus_dic.values())
            min_nodes = 3

            for cluster in set(clusters):
                if counts[cluster] > min_nodes:
                    att_frequency_clus[cluster] = {}
                    for code in self.codes:
                        att_value = (self.node_attributes[code][attribute])
                        if clus_dic[code] == cluster:
                            if att_value not in att_frequency_clus[cluster]:
                                att_frequency_clus[cluster][att_value]  = 1
                            else:
                                att_frequency_clus[cluster][att_value]  += 1
                        elif att_value not in  att_frequency_clus[cluster]:
                            att_frequency_clus[cluster][att_value]  = 0


            def sdi(data):
                def p(n, N):
                    if n == 0:
                        return 0
                    else:
                        return ((float(n)/N) * ln(float(n)/N))

                N = sum(data.values())

                return -sum(p(n, N) for n in data.values() if n != 0)



            shannon = {}
            for cluster,att in att_frequency_clus.items():
                shannon[cluster] = sdi(att)

            print(shannon)

            results[threshold] = (sum(list(shannon.values())))

        return(results)


    def plotSimilarityMatrix(self):
        """
        Plots the similarity matrix of the given sequences
        """

        plt.matshow(self.similarity_matrix, vmin=0.0, vmax=1.0)
        cbar = plt.colorbar()
        cbar.set_label('Sequence identity')
        plt.xlabel('Sequence index')
        plt.ylabel('Sequence index')
