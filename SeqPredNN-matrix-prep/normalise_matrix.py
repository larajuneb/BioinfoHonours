from Bio.PDB import *
import os
import sys

path = sys.argv[1]
mode = ""
if "oversample" in path:
    mode = "oversampled"
if "undersample" in path:
    mode = "undersampled"

unnormalised_matrix = [[0 for x in range(20)] for y in range(20)] #[[true HIS, (true res index),(pred res index)][true ARG][true LYS][etc ...]]
normalised_matrix = [[0 for x in range(20)] for y in range(20)] #[[true HIS, (true res index),(pred res index)][true ARG][true LYS][etc ...]]
rounded_normalised_matrix = [[0 for x in range(20)] for y in range(20)] #[[true HIS, (true res index),(pred res index)][true ARG][true LYS][etc ...]]

predicted_res = ""
true_res = ""
res_dict = {"HIS": 0, "ARG": 1, "LYS": 2, "GLN": 3, "GLU": 4, "ASP": 5, "ASN": 6, "GLY": 7, "ALA": 8, "SER": 9, "THR": 10, "PRO": 11, "CYS": 12, "VAL": 13, "ILE": 14, "MET": 15, "LEU": 16, "PHE": 17, "TYR": 18, "TRP": 19}
support_dict = {"HIS": 0, "ARG": 0, "LYS": 0, "GLN": 0, "GLU": 0, "ASP": 0, "ASN": 0, "GLY": 0, "ALA": 0, "SER": 0, "THR": 0, "PRO": 0, "CYS": 0, "VAL": 0, "ILE": 0, "MET": 0, "LEU": 0, "PHE": 0, "TYR": 0, "TRP": 0}

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

# build normalised matrix
for row in range(20):
    res_name = list(res_dict)[row]
    line = res_name + ": "
    total_check = 0
    for col in range (20):
        normalised_matrix[row][col] = unnormalised_matrix[row][col] / support_dict.get(res_name)
        rounded_normalised_matrix[row][col] = round(unnormalised_matrix[row][col] / support_dict.get(res_name), 4)
        line += str(round(normalised_matrix[row][col], 3)) + ", "
        total_check += normalised_matrix[row][col]

# build csv
filename = mode + "_normalised_matrix.csv"
with open(filename, mode="w") as file:
    file.write(" ,HIS,ARG,LYS,GLN,GLU,ASP,ASN,GLY,ALA,SER,THR,PRO,CYS,VAL,ILE,MET,LEU,PHE,TYR,TRP\n") #first line
    counter = 0
    for res_list in normalised_matrix:
        file.write(list(res_dict)[counter] + "," + ",".join(map(str, res_list)) + "\n")
        counter += 1

filename = mode + "_rounded_normalised_matrix.csv"
with open(filename, mode="w") as file:
    file.write(" ,HIS,ARG,LYS,GLN,GLU,ASP,ASN,GLY,ALA,SER,THR,PRO,CYS,VAL,ILE,MET,LEU,PHE,TYR,TRP\n") #first line
    counter = 0
    for res_list in rounded_normalised_matrix:
        file.write(list(res_dict)[counter] + " (true)," + ",".join(map(str, res_list)) + "\n")
        counter += 1