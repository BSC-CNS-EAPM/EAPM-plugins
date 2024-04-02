import argparse
import json
import MDAnalysis
import scipy
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--protein")
parser.add_argument("-t", "--threshold", default=3.5)

args = parser.parse_args()
protein = args.protein
threshold = float(args.threshold)

def search_catalytic_triad(prot, threshold=3.5):

	serine_og = prot.select_atoms("(resname SER and name OG)")
	histidine_nx = prot.select_atoms("(resname HIS and (name NE2 or name ND1))")
	acid_ox = prot.select_atoms("(resname ASP and name OD1 or name OD2) or (resname GLU and name OE1 or name OE2)")

	coors1 = serine_og.positions
	coors2 = histidine_nx.positions
	coors3 = acid_ox.positions

	a = scipy.spatial.distance_matrix(coors2, coors1)
	b = scipy.spatial.distance_matrix(coors2, coors3)

	a_clos = (np.argwhere(a<threshold))
	b_clos = (np.argwhere(b<threshold))

	dict_a_clos = {}

	for x,y in a_clos:
	    if x not in dict_a_clos:
	        dict_a_clos[x] = [y]
	    else:
	        dict_a_clos[x].append(y)

	dict_b_clos = {}

	for x,y in b_clos:
	    if x not in dict_b_clos:
	        dict_b_clos[x] = [y]
	    else:
	        dict_b_clos[x].append(y)

	ser = []
	his = []
	asp = []

	for his_atom1,ser_list in dict_a_clos.items():
		for his_atom2,asp_list in dict_b_clos.items():
			his1 = histidine_nx[his_atom1]
			his2 = histidine_nx[his_atom2]
			if his1.resid == his2.resid and his1.chainID == his2.chainID and his1 != his2:
				for ser_atom in ser_list:
					for asp_atom in asp_list:
						ser.append(ser_atom)
						asp.append(asp_atom)
						his.append(his_atom1)


	ser_names = [(serine_og[x].chainID,serine_og[x].resname,str(serine_og[x].resid)) for x in (ser)]
	his_names = [(histidine_nx[x].chainID,histidine_nx[x].resname,str(histidine_nx[x].resid),histidine_nx[x].name) for x in (his)]
	asp_names = [(acid_ox[x].chainID,acid_ox[x].resname,str(acid_ox[x].resid)) for x in (asp)]

	res_triads = [x for x in zip(ser_names, his_names, asp_names)]
	res_triads = {i:x for i,x in enumerate(list(set(res_triads)))}


	return(res_triads)

prot = MDAnalysis.Universe(protein)
triads = (search_catalytic_triad(prot, threshold=threshold))

if triads != {}:
	with open(protein.split("/")[-1] + ".json","w") as jf:
		json.dump(triads, jf)
else:
	print("Protein with code "+str(protein)+" has no triad")
