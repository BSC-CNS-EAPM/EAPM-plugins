import os
import io

import json
import argparse

import gzip
from Bio.PDB import PDBParser

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--protein")
args = parser.parse_args()
protein = args.protein


def pdb_chain_dictionary(protein):
    chain_ids = {}
    file = gzip.open(protein,"rt")
    struct = PDBParser().get_structure("protein", file)
    chains = struct[0]
    for i,chain in enumerate(chains):
        chain_ids[i] = (chain.get_id())
    return (chain_ids)


with open(protein.split("/")[-1] + ".json","w") as jf:
    json.dump(pdb_chain_dictionary(protein), jf)
