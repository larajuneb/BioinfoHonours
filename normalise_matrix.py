from Bio.PDB import *
import os
import sys

path = sys.argv[1]

unnormalised_matrix = [[0 for x in range(20)] for y in range(20)] #[[true HIS, (true res index),(pred res index)][true ARG][true LYS][etc ...]]
normalised_matrix = [[0 for x in range(20)] for y in range(20)] #[[true HIS, (true res index),(pred res index)][true ARG][true LYS][etc ...]]

predicted_res = ""
true_res = ""
res_dict = {"HIS": 0, "ARG": 1, "LYS": 2, "GLN": 3, "GLU": 4, "ASP": 5, "ASN": 6, "GLY": 7, "ALA": 8, "SER": 9, "THR": 10, "PRO": 11, "CYS": 12, "VAL": 13, "ILE": 14, "MET": 15, "LEU": 16, "PHE": 17, "TYR": 18, "TRP": 19}
support_dict = {"HIS": 0, "ARG": 1, "LYS": 2, "GLN": 3, "GLU": 4, "ASP": 5, "ASN": 6, "GLY": 7, "ALA": 8, "SER": 9, "THR": 10, "PRO": 11, "CYS": 12, "VAL": 13, "ILE": 14, "MET": 15, "LEU": 16, "PHE": 17, "TYR": 18, "TRP": 19}

# read in raw data from csv
with open(path) as file:
    for line in file:
        row = line.rstrip()
        if not row.startswith("Predicted"):
            predicted_res = row[:3]
            predicted_index = res_dict.get(predicted_res)
            true_res = row[4:7]
            true_index = res_dict.get(true_res)
            value = int(row[8:])
            unnormalised_matrix[predicted_index][true_index] = value

# build unnormalised matrix
print("UNNORMALISED MATRIX: ")
line = ""
for row in range(20):
    line = ""
    line = list(res_dict)[row] + ": "
    res_support = 0
    for col in range (20):
        line += str(unnormalised_matrix[row][col]) + ", "
        res_support += unnormalised_matrix[row][col]
    support_dict[list(res_dict)[row]] = res_support
    line += "==(" + str(res_support) + ")"
    print(line)

print()
print(support_dict)
print()

# build normalised matrix
print("NORMALISED MATRIX: ")
for row in range(20):
    res_name = list(res_dict)[row]
    line = res_name + ": "
    total_check = 0
    for col in range (20):
        normalised_matrix[row][col] = unnormalised_matrix[row][col] / support_dict.get(res_name)
        line += str(round(normalised_matrix[row][col], 3)) + ", "
        total_check += normalised_matrix[row][col]
    line += "==(" + str(total_check) + ")"
    print(line)

# build csv, example code: 
# with open("new_SeqPredNN_pdb_subset.csv", mode="w") as file:
#     file.write("Protein,Filename,Chain\n") #first line
#     for entry in official_identifiers_list:
#         file.write(str(entry[0]) + "," + str(entry[0]) + ".pdb.gz," + str(entry[1]) + "\n")