from Bio.PDB import *
import os

pdb_list = PDBList()
ids = []
chains = []
pdb_dict = {}

with open("/home/larajuneb/Honours/PROJECT_(721)/Coding/SeqPredNN/SeqPredNN-main/Hons-Project/cullpdb_pc90.0_res0.0-2.2_len40-10000_R0.25_Xray_d2023_07_18_chains31896") as file:
    for line in file:
        row = line.rstrip()
        if not row.startswith("PDBchain"):
            identifier = row[:row.find(" ")]
            id = identifier[:4]
            chain = identifier[4:]
            ids.append(id)
            chains.append(chain)
            pdb_dict[id] = chain

directory = "/home/larajuneb/Honours/PROJECT_(721)/Coding/SeqPredNN/SeqPredNN-main/Hons-Project/new_protein_data_for_training/PDB_files"
unique_ids = set(ids)

# print(len(ids))

# loop through files in directory, if contained in IDs list, remove from IDs list
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        if filename[3:7] in unique_ids:
            # print("yes")
            # print(filename[3:7])
            unique_ids.remove(filename[3:7])

# print(len(ids))


pdb_list.download_pdb_files(unique_ids, pdir=directory, file_format="pdb")