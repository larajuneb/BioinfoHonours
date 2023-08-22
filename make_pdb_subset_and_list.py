import os
import random

identifiers = []
ids_to_download = []

# make list of all PDB IDs to download
with open("/home/larajuneb/Honours/PROJECT_(721)/Coding/SeqPredNN/SeqPredNN-main/Hons-Project/cullpdb_pc90.0_res0.0-2.2_len40-10000_R0.25_Xray_d2023_07_18_chains31896") as file:
    for line in file:
        row = line.rstrip()
        if not row.startswith("PDBchain"):
            identifier = row[:row.find(" ")]
            id = identifier[:4]
            chain = identifier[4:]
            ids_to_download.append(id)
            identifiers.append([id, chain])

print("Total length of identifiers list: " + str(len(identifiers)))

# make a list of all PDB IDs successfully downloaded
ids_downloaded = []
directory = "/home/larajuneb/Honours/PROJECT_(721)/Coding/SeqPredNN/SeqPredNN-main/Hons-Project/new_protein_data_for_training/PDB_files_batch"
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        ids_downloaded.append(filename[:4])

print("Number of PDB files actually downloaded: " + str(len(ids_downloaded)))

official_identifiers_list = []
official_id_list_check = []

for pair in identifiers:
    if pair[0] in ids_downloaded:
        official_identifiers_list.append(pair)
        official_id_list_check.append(pair[0])

undownloaded_identifiers = []

# make a list of all idenitifiers not downloaded (ID and chain)
for set in identifiers:
    if not set[0] in ids_downloaded:
        undownloaded_identifiers.append(set)

print("Length of official identifiers after removing non-downloaded PDB IDs: " + str(len(official_identifiers_list)))
print("Length of identifiers of non-downloaded PDB IDs: " + str(len(undownloaded_identifiers)))

if (len(undownloaded_identifiers) + len(official_identifiers_list)) == len(identifiers):
    print("TRUE: " + str(len(undownloaded_identifiers) + len(official_identifiers_list)))

with open("new_SeqPredNN_pdb_subset.csv", mode="w") as file:
    file.write("Protein,Filename,Chain\n") #first line
    for entry in official_identifiers_list:
        file.write(str(entry[0]) + "," + str(entry[0]) + ".pdb.gz," + str(entry[1]) + "\n")

official_chains = []
for entry in official_identifiers_list:
    official_chains.append(str(entry[0])+str(entry[1]))

n = round(0.1*len(official_chains))
print(n)
test_set_chains = random.sample(official_chains, n)

# make random list of 10% chains for test set
with open("test_set.txt", mode="w") as file:
    for test_chain in test_set_chains:
        file.write(str(test_chain) + "\n")