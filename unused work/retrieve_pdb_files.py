from Bio.PDB import *
# import time
# import os

# def download_pdb(pdb_id: str, directory: str) -> None:
#     return

pdb_list = PDBList()
ids = []
chains = []
pdb_dict = {}

with open("/home/larajuneb/Honours/PROJECT_(721)/Coding/SeqPredNN/SeqPredNN-main/Hons-Project/new_protein_data_for_training/cullpdb_pc90.0_res0.0-2.2_len40-10000_R0.25_Xray_d2023_07_18_chains31896") as file:
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
# print("Total number of IDs: " + str(len(ids)))
# print("Total number of unique IDs: " + str(len(unique_ids)))
# print("Length of dictionary: " + str(len(pdb_dict)))

pdb_list.download_pdb_files(unique_ids, pdir=directory, file_format="pdb")

# start = time.time()
# for filename in os.listdir(directory):
#     f = os.path.join(directory, filename)
#     # checking if it is a file
#     if os.path.isfile(f):
#         print(filename)
# end = time.time()
# print(str(((end-start)/10)*31896))

# download all files
# if number of files is equal to number of keys in the dictionary, then all files were downloaded
# and you can add all the info from the list to the csv