import os
import io
import bioprospecting
import prepare_proteins
import json
import argparse
import shutil
import random
from Bio import SeqIO
from Bio import PDB
import gzip


pdb_list = []
#for filename in sorted(os.listdir("catalytic_triads"))[start:end]:
for filename in os.listdir("catalytic_triads"):
    if filename.endswith(".json"):
        pdb_list.append(filename[3:7])

pdb_data = bioprospecting.databases.scrape.getPDBInformation(pdb_list,'triad_classification/pdb_data_')

up_pdb_dic = {}
pdb_code_dic = {}
up_data = {}

for code in pdb_list:
    if code not in pdb_code_dic:
        pdb_code_dic[code] = {}
    for chain in pdb_data[code]['Macromolecules']:
        up_code = pdb_data[code]['Macromolecules'][chain]['Uniprot']

        if len(chain) != 1:
            chain = chain[-2]

        if up_code is not None:
            if up_code not in up_pdb_dic:
                up_pdb_dic[up_code] = []
            up_pdb_dic[up_code].append(code+'_'+chain)

            if os.path.exists("triad_classification/up_data_"+up_code+".json"):
                output_file = open("triad_classification/up_data_"+up_code+".json")
                single_up_data = json.load(output_file)
                up_data[up_code] = single_up_data
                output_file.close()
                print('UniProt data found for code %s.' % str(up_code))

            else:
                up_data[up_code] = bioprospecting.databases.getUniprotData(up_code)
                with open("triad_classification/up_data_"+up_code+".json","w") as jf:
                    print('UniProt data file for %s written.' % str(up_code))
                    json.dump(up_data[up_code], jf)

        pdb_code_dic[code][chain] = up_code

with open("up_pdb_dic.json","w") as jf:
    json.dump(up_pdb_dic, jf)

real_triads = {}
dif_real_triads = {}

doc_triads = {}
found_triads = {}

ambigous_triads = {}

for filename in os.listdir("catalytic_triads"):
#for filename in sorted(os.listdir("catalytic_triads"))[start:end]:
    if filename.endswith(".json"):
        #doc_triads = []
        #found_triads = []
        with open('catalytic_triads/'+filename, 'r') as f:
            for index,found_triad in json.load(f).items():
                real_triad_check = False
                pos_found_triad = set((found_triad[0][2],found_triad[1][2],found_triad[2][2]))
                #chain = chain_name_dic[filename][str(found_triad[0][0])]
                chain = str(found_triad[0][0])
                if found_triad[0][0] == found_triad[1][0] == found_triad[2][0]:
                    if chain in pdb_code_dic[filename[3:7]]:
                        if pdb_code_dic[filename[3:7]][chain] is not None:
                            if up_data[pdb_code_dic[filename[3:7]][chain]] != None:
                                if "active site" in up_data[pdb_code_dic[filename[3:7]][chain]]["features"]:
                                    real_triad = set(up_data[pdb_code_dic[filename[3:7]][chain]]["features"]["active site"])

                                    ## Get pdb positions of uniprot triad
                                    pdb_positions = []
                                    with gzip.open('/home/bscuser/Databases/PDB/'+filename[:-5],'rt') as pdb_file:
                                        #for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                                        parser = PDB.PDBParser()
                                        structure = parser.get_structure("model",pdb_file)
                                        model = structure[0]
                                        for pdb_chain in model:
                                            if chain == pdb_chain.id and 'Sequence' in up_data[pdb_code_dic[filename[3:7]][chain]]:
                                                up_seq = up_data[pdb_code_dic[filename[3:7]][chain]]['Sequence']
                                                pdb_seq = ''
                                                residue_ids = {'pdb':{}}
                                                for i,r in enumerate(pdb_chain.get_residues()):
                                                    if r.id[0] == ' ': # Non heteroatom filter
                                                        try:
                                                            pdb_seq += PDB.Polypeptide.three_to_one(r.resname)
                                                        except:
                                                            pdb_seq += 'X'
                                                    else:
                                                        pdb_seq += 'X'

                                                    residue_ids['pdb'][i+1] = r.id[1]

                                                print(filename,chain)
                                                print(pdb_seq)
                                                print('#####')
                                                print(up_seq)

                                                to_align = {'pdb':pdb_seq,'up':up_seq}
                                                aln = prepare_proteins.alignment.mafft.multipleSequenceAlignment(to_align)

                                                residue_positions = {}

                                                residue_positions['pdb'] = 0
                                                residue_positions['up'] = 0

                                                ## Check chain names are consistent for PDB UniProt

                                                for i in range(aln.get_alignment_length()):
                                                    for entry in aln:
                                                        if entry.seq[i] != '-':
                                                            residue_positions[entry.id] += 1

                                                    if str(residue_positions['up']) in real_triad and residue_positions['pdb'] in residue_ids['pdb']:
                                                        pdb_positions.append(str(residue_ids['pdb'][residue_positions['pdb']]))

                                    print('up',real_triad)
                                    print('pdb',pos_found_triad)
                                    print('new_positions',pdb_positions)

                                    if pdb_positions != []:
                                        real_triad = set(pdb_positions)

                                    if pos_found_triad == real_triad:
                                        real_triad_check = True
                                        #doc_triads.append(found_triad)
                                        if filename not in doc_triads:
                                            doc_triads[filename] = []
                                        doc_triads[filename].append(index)

                                    elif real_triad.issubset(pos_found_triad) or real_triad.issuperset(pos_found_triad):
                                        real_triad_check = True
                                        #if filename not in doc_triads:
                                        #    doc_triads[filename] = []
                                        #doc_triads[filename].append(index)

                                        if filename not in ambigous_triads:
                                            ambigous_triads[filename] = []
                                        ambigous_triads[filename].append(index)

                                        #if filename not in real_triads:
                                        #    real_triads[filename] = []
                                        #real_triads[filename][chain] = (pdb_code_dic[filename[3:7]][chain])
                                        #real_triads[filename].append(found_triad)

                                    #else:
                                    #    real_triad_check = True
                                    #    if filename not in dif_real_triads:
                                    #        dif_real_triads[filename] = []
                                    #    dif_real_triads[filename].append(found_triad)

                if real_triad_check == False and (found_triad[0][0] == found_triad[1][0] == found_triad[2][0]):
                    #found_triads.append(found_triad)
                    if filename not in found_triads:
                        found_triads[filename] = []
                    found_triads[filename].append(index)

                    #if filename not in found_triads:
                    #    found_triads[filename] = []
                    #found_triads[filename][chain] = (pdb_code_dic[filename[3:7]][chain])
                    #found_triads[filename].append(found_triad)

        # if not os.path.exists("doc_triads"):
        #     os.mkdir("doc_triads")
        #
        # if not os.path.exists("found_triads"):
        #     os.mkdir("found_triads")
        #
        # if doc_triads != []:
        #     with open("doc_triads/"+filename,"w") as jf:
        #         json.dump(doc_triads, jf)
        #
        # if found_triads != []:
        #     with open("found_triads/"+filename,"w") as jf:
        #         json.dump(found_triads, jf)





#num_doc = len(next(os.walk("doc_triads/"))[2])

#for filename in random.sample(os.listdir("found_triads/"),num_doc):
#    shutil.move(os.path.join("found_triads/",filename),"sample_found_triads/")

with open("doc_triads.json","w") as jf:
    json.dump(doc_triads, jf)

with open("found_triads.json","w") as jf:
    json.dump(found_triads, jf)

with open("ambigous_triads.json","w") as jf:
    json.dump(ambigous_triads, jf)

#with open("error_doc_triads.json","w") as jf:
#    json.dump(dif_real_triads, jf)
