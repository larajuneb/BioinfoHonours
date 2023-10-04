from Bio.PDB import *
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

path = sys.argv[1]
mode = ""
if "oversample" in path:
    mode = "oversampled"
if "undersample" in path:
    mode = "undersampled"

if len(sys.argv) < 3:
    print("Please input output mode (0, ..., 7) as 2nd argument.")
    exit(0)

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

# codons_per_aa = {"HIS": ['CAU', 'CAC'], "ARG": ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], "LYS": ['AAA', 'AAG'], "GLN": ['CAA', 'CAG'], "GLU": ['GAA', 'GAG'], "ASP": ['GAU', 'GAC'], "ASN": ['AAU', 'AAC'], "GLY": ['GGU', 'GGC', 'GGA', 'GGG'], "ALA": ['GCU', 'GCC', 'GCA', 'GCG'], "SER": ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], "THR": ['ACU', 'ACC', 'ACA', 'ACG'], "PRO": ['CCU', 'CCC', 'CCA', 'CCG'], "CYS": ['UGU', 'UGC'], "VAL": ['GUU', 'GUC', 'GUA', 'GUG'], "ILE": ['AUU', 'AUC', 'AUA'], "MET": ['AUG'], "LEU": ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'], "PHE": ['UUU', 'UUC'], "TYR": ['UAU', 'UAC'], "TRP": ['UGG']}

codons_per_aa = {"HIS": [['C'], ['A'], ['U', 'C']], "ARG": [[['C'], ['G'], ['U', 'C', 'A', 'G']],[['A'], ['G'], ['A', 'G']]], "LYS": [['A'], ['A'], ['A', 'G']], "GLN": [['C'], ['A'], ['A', 'G']], "GLU": [['G'], ['A'], ['A', 'G']], "ASP": [['G'], ['A'], ['U', 'C']], "ASN": [['A'], ['A'], ['U', 'C']], "GLY": [['G'], ['G'], ['U', 'C', 'A', 'G']], "ALA": [['G'], ['C'], ['U', 'C', 'A', 'G']], "SER": [[['U'], ['C'], ['U', 'C', 'A', 'G']], [['A'], ['G'], ['U', 'C']]], "THR": [['A'], ['C'], ['U', 'C', 'A', 'G']], "PRO": [['C'], ['C'], ['U', 'C', 'A', 'G']], "CYS": [['U'], ['G'], ['U', 'C']], "VAL": [['G'], ['U'], ['U', 'C', 'A', 'G']], "ILE": [['A'], ['U'], ['U', 'C', 'A']], "MET": [['A'], ['U'], ['G']], "LEU": [[['C'], ['U'], ['U', 'C', 'A', 'G']], [['U'], ['U'], ['A', 'G']]], "PHE": [['U'], ['U'], ['U', 'C']], "TYR": [['U'], ['A'], ['U', 'C']], "TRP": [['U'], ['G'], ['G']]}

codon_comparison_table = [['0' for x in range(20)] for y in range(20)] #number of nucleotide changes necessary, y axis (row) = true, x axis (col) = pred
codon_position_changes = {1:0, 2:0, 3:0} #increase to 1 if change in that position is necessary
position_matrix = [['0' for x in range(20)] for y in range(20)] #position in codon of mutations necessary

treat_different = ['LEU', 'ARG', 'SER']

for true in range(20): #row
    for pred in range(20): #col
        #if name is LEU, ARG, or SER, look at 2 different sets of codons
        true_name = list(res_dict)[true]
        pred_name = list(res_dict)[pred]
        true_codons = codons_per_aa.get(true_name)
        pred_codons = codons_per_aa.get(pred_name)
        positions = []
        
        if true_name not in treat_different and pred_name not in treat_different:
            true_first = true_codons[0]
            true_second = true_codons[1]
            true_third = true_codons[2]

            pred_first = pred_codons[0]
            pred_second = pred_codons[1]
            pred_third = pred_codons[2]

            diff1 = list(set(true_first).symmetric_difference(set(pred_first)))
            diff2 = list(set(true_second).symmetric_difference(set(pred_second)))
            diff3 = list(set(true_third).symmetric_difference(set(pred_third)))

            if len(diff1) != 0:
                codon_position_changes[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                codon_position_changes[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                codon_position_changes[3] = 1
                positions.append(3)

            codon_comparison_table[true][pred] = str(sum(list(codon_position_changes.values())))

        if true_name in treat_different and pred_name not in treat_different:
            set1_diff = {1:0, 2:0, 3:0}
            set2_diff = {1:0, 2:0, 3:0}
            set1 = true_codons[0]
            set2 = true_codons[1]

            pred_first = pred_codons[0]
            pred_second = pred_codons[1]
            pred_third = pred_codons[2]

            diff1 = list(set(set1[0]).symmetric_difference(set(pred_first)))
            diff2 = list(set(set1[1]).symmetric_difference(set(pred_second)))
            diff3 = list(set(set1[2]).symmetric_difference(set(pred_third)))

            if len(diff1) != 0:
                set1_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set1_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set1_diff[3] = 1
                positions.append(3)

            diff1 = list(set(set2[0]).symmetric_difference(set(pred_first)))
            diff2 = list(set(set2[1]).symmetric_difference(set(pred_second)))
            diff3 = list(set(set2[2]).symmetric_difference(set(pred_third)))

            if len(diff1) != 0:
                set2_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set2_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set2_diff[3] = 1
                positions.append(3)
            
            diffs = [sum(list(set1_diff.values())), sum(list(set2_diff.values()))]
            diffs_np = np.array(diffs)
            diffs_set = np.unique(diffs_np)
            diffs_set.sort()

            if len(diffs_set) > 1:
                codon_comparison_table[true][pred] = str(diffs_set[0]) + "|" + str(diffs_set[1])
            else:
                codon_comparison_table[true][pred] = str(diffs_set[0])

        if true_name not in treat_different and pred_name in treat_different:
            set1_diff = {1:0, 2:0, 3:0}
            set2_diff = {1:0, 2:0, 3:0}
            set1 = pred_codons[0]
            set2 = pred_codons[1]

            true_first = true_codons[0]
            true_second = true_codons[1]
            true_third = true_codons[2]

            diff1 = list(set(set1[0]).symmetric_difference(set(true_first)))
            diff2 = list(set(set1[1]).symmetric_difference(set(true_second)))
            diff3 = list(set(set1[2]).symmetric_difference(set(true_third)))

            if len(diff1) != 0:
                set1_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set1_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set1_diff[3] = 1
                positions.append(3)

            diff1 = list(set(set2[0]).symmetric_difference(set(true_first)))
            diff2 = list(set(set2[1]).symmetric_difference(set(true_second)))
            diff3 = list(set(set2[2]).symmetric_difference(set(true_third)))

            if len(diff1) != 0:
                set2_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set2_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set2_diff[3] = 1
                positions.append(3)
            
            diffs = [sum(list(set1_diff.values())), sum(list(set2_diff.values()))]
            diffs_np = np.array(diffs)
            diffs_set = np.unique(diffs_np)
            diffs_set.sort()

            if len(diffs_set) > 1:
                codon_comparison_table[true][pred] = str(diffs_set[0]) + "|" + str(diffs_set[1])
            else:
                codon_comparison_table[true][pred] = str(diffs_set[0])

        if true_name in treat_different and pred_name in treat_different:
            #do something
            # a b
            # c d
            set13_diff = {1:0, 2:0, 3:0} #ac
            set14_diff = {1:0, 2:0, 3:0} #ad
            set23_diff = {1:0, 2:0, 3:0} #bc
            set24_diff = {1:0, 2:0, 3:0} #bd
            set1 = true_codons[0] #a
            set2 = true_codons[1] #b
            set3 = pred_codons[0] #c
            set4 = pred_codons[1] #d

            # SET 13
            diff1 = list(set(set1[0]).symmetric_difference(set(set3[0])))
            diff2 = list(set(set1[1]).symmetric_difference(set(set3[1])))
            diff3 = list(set(set1[2]).symmetric_difference(set(set3[2])))

            if len(diff1) != 0:
                set13_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set13_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set13_diff[3] = 1
                positions.append(3)

            # SET 14
            diff1 = list(set(set1[0]).symmetric_difference(set(set4[0])))
            diff2 = list(set(set1[1]).symmetric_difference(set(set4[1])))
            diff3 = list(set(set1[2]).symmetric_difference(set(set4[2])))

            if len(diff1) != 0:
                set14_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set14_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set14_diff[3] = 1
                positions.append(3)

            # SET 23
            diff1 = list(set(set2[0]).symmetric_difference(set(set3[0])))
            diff2 = list(set(set2[1]).symmetric_difference(set(set3[1])))
            diff3 = list(set(set2[2]).symmetric_difference(set(set3[2])))

            if len(diff1) != 0:
                set23_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set23_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set23_diff[3] = 1
                positions.append(2)

            # SET 24
            diff1 = list(set(set2[0]).symmetric_difference(set(set4[0])))
            diff2 = list(set(set2[1]).symmetric_difference(set(set4[1])))
            diff3 = list(set(set2[2]).symmetric_difference(set(set4[2])))

            if len(diff1) != 0:
                set24_diff[1] = 1
                positions.append(1)
            if len(diff2) != 0:
                set24_diff[2] = 1
                positions.append(2)
            if len(diff3) != 0:
                set24_diff[3] = 1
                positions.append(3)
            
            diffs = [sum(list(set13_diff.values())), sum(list(set14_diff.values())), sum(list(set23_diff.values())), sum(list(set24_diff.values()))]
            diffs_np = np.array(diffs)
            diffs_set = np.unique(diffs_np)
            diffs_set.sort()

            if len(diffs_set) > 1:
                if diffs_set[0] == 0:
                    codon_comparison_table[true][pred] = str(diffs_set[1])
                else:
                    codon_comparison_table[true][pred] = str(diffs_set[0]) + "|" + str(diffs_set[1])
            else:
                codon_comparison_table[true][pred] = str(diffs_set[0])
            
            if true_name == pred_name:
                codon_comparison_table[true][pred] = '0'

        positions_np = np.array(positions)
        positions_set = np.unique(positions_np)
        positions_set.sort()
        if len(positions) == 0 or true_name == pred_name:
            position_matrix[true][pred] = "0"
        else:
            position_matrix[true][pred] = "|".join(map(str, positions_set))
        codon_position_changes[1] = 0
        codon_position_changes[2] = 0
        codon_position_changes[3] = 0

if sys.argv[2] == '0': # build csv
    filename = "codon_comparison_matrix.csv"
    with open(filename, mode="w") as file:
        file.write(" ,HIS,ARG,LYS,GLN,GLU,ASP,ASN,GLY,ALA,SER,THR,PRO,CYS,VAL,ILE,MET,LEU,PHE,TYR,TRP\n") #first line
        counter = 0
        for res_list in codon_comparison_table:
            file.write(list(res_dict)[counter] + "," + ",".join(map(str, res_list)) + "\n")
            counter += 1

    filename = "positions_matrix.csv"
    with open(filename, mode="w") as file:
        file.write(" ,HIS,ARG,LYS,GLN,GLU,ASP,ASN,GLY,ALA,SER,THR,PRO,CYS,VAL,ILE,MET,LEU,PHE,TYR,TRP\n") #first line
        counter = 0
        for res_list in position_matrix:
            file.write(list(res_dict)[counter] + "," + ",".join(map(str, res_list)) + "\n")
            counter += 1

#creating plots

amino_acids = ["H", "R", "K", "Q", "E", "D", "N", "G", "A", "S", "T", "P", "C", "V", "I", "M", "L", "F", "Y", "W"]

x = []
y = []
z = []
dx = []
dy = []
dz = []
colours = []

# ---------------------------------------------1) FULL BAR PLOT, Z AXIS = POSITION OF CHANGE, COLOURS = # OF CHANGES---------------------------------------------------------------
if sys.argv[2] == '1':
    for row in range(20):
        for col in range(20):
            x.append(row)
            y.append(col)
            z.append(0)
            dx.append(1)
            dy.append(1)
            if '|' not in position_matrix[row][col]:
                dz.append(int(position_matrix[row][col]))
            else:
                if len(position_matrix[row][col]) == 5:
                    dz.append((int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5]))/3)
                elif len(position_matrix[row][col]) == 3:
                    dz.append((int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]))/2)

            if '|' in codon_comparison_table[row][col]:
                i = (int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3])) / 2
                if i == 1.5:
                    colours.append("turquoise")
                elif i == 2:
                    colours.append("red") #1/3
                elif i == 2.5:
                    colours.append("yellow") #2/3

            elif int(codon_comparison_table[row][col]) == 0:
                colours.append("magenta")
            elif int(codon_comparison_table[row][col]) == 1:
                colours.append("blue")
            elif int(codon_comparison_table[row][col]) == 2:
                colours.append("green")
            elif int(codon_comparison_table[row][col]) == 3:
                colours.append("orange")

    center = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]


    fig = plt.figure(figsize=(15, 9))
    ax = plt.axes(projection="3d")
    ax.bar3d(x, y, z, dx, dy, dz, color=colours)
    ax.set_xlabel("True amino acid")
    ax.set_ylabel("Predicted amino acid")
    ax.set_zlabel("Position of change")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("FULL BAR PLOT, Z AXIS = POSITION OF CHANGE, COLOURS = # OF CHANGES")

    magenta_patch = mpatches.Patch(color='magenta', label='# changes = 0')
    blue_patch = mpatches.Patch(color='blue', label='# changes = 1')
    turq_patch = mpatches.Patch(color='turquoise', label='# changes = 1/2')
    red_patch = mpatches.Patch(color='red', label='# changes = 1/3')
    green_patch = mpatches.Patch(color='green', label='# changes = 2')
    yellow_patch = mpatches.Patch(color='yellow', label='# changes = 2/3')
    orange_patch = mpatches.Patch(color='orange', label='# changes = 3')
    plt.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch], loc='upper left', framealpha=0.5, title="Position of change")

    plt.show()
# ---------------------------------------------END---------------------------------------------------------------

# ---------------------------------------------2) HALF BAR PLOT, Z AXIS = POSITION OF CHANGE, COLOURS = # OF CHANGES---------------------------------------------------------------
if sys.argv[2] == '2':
    for row in range(20):
        for col in range(row, 20):
            x.append(row)
            y.append(col)
            z.append(0)
            dx.append(1)
            dy.append(1)
            if '|' not in position_matrix[row][col]:
                dz.append(int(position_matrix[row][col]))
            else:
                if len(position_matrix[row][col]) == 5:
                    dz.append((int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5]))/3)
                elif len(position_matrix[row][col]) == 3:
                    dz.append((int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]))/2)

            if '|' in codon_comparison_table[row][col]:
                i = (int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3])) / 2
                if i == 1.5:
                    colours.append("turquoise")
                elif i == 2:
                    colours.append("red") #1/3
                elif i == 2.5:
                    colours.append("yellow") #2/3

            elif int(codon_comparison_table[row][col]) == 0:
                colours.append("magenta")
            elif int(codon_comparison_table[row][col]) == 1:
                colours.append("blue")
            elif int(codon_comparison_table[row][col]) == 2:
                colours.append("green")
            elif int(codon_comparison_table[row][col]) == 3:
                colours.append("orange")

    center = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]


    fig = plt.figure(figsize=(15, 9))
    ax = plt.axes(projection="3d")
    ax.bar3d(x, y, z, dx, dy, dz, color=colours)
    ax.set_xlabel("True amino acid")
    ax.set_ylabel("Predicted amino acid")
    ax.set_zlabel("Position of change")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("HALF BAR PLOT, Z AXIS = POSITION OF CHANGE, COLOURS = # OF CHANGES")

    magenta_patch = mpatches.Patch(color='magenta', label='# changes = 0')
    blue_patch = mpatches.Patch(color='blue', label='# changes = 1')
    turq_patch = mpatches.Patch(color='turquoise', label='# changes = 1/2')
    red_patch = mpatches.Patch(color='red', label='# changes = 1/3')
    green_patch = mpatches.Patch(color='green', label='# changes = 2')
    yellow_patch = mpatches.Patch(color='yellow', label='# changes = 2/3')
    orange_patch = mpatches.Patch(color='orange', label='# changes = 3')
    plt.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch], loc='upper left', framealpha=0.5, title="Number of changes")

    plt.show()

# ---------------------------------------------END---------------------------------------------------------------

# ---------------------------------------------3) FULL BAR PLOT, Z AXIS = # OF CHANGES, COLOURS = POSITION OF CHANGE---------------------------------------------------------------
if sys.argv[2] == '3':
    for row in range(20):
        for col in range(20):
            x.append(row)
            y.append(col)
            z.append(0)
            dx.append(1)
            dy.append(1)
            if '|' not in codon_comparison_table[row][col]:
                dz.append(int(codon_comparison_table[row][col]))
            else:
                if len(codon_comparison_table[row][col]) == 5:
                    dz.append((int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]) + int(codon_comparison_table[row][col][4:5]))/3)
                elif len(codon_comparison_table[row][col]) == 3:
                    dz.append((int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]))/2)

            if '|' in position_matrix[row][col]:
                i = 0
                if len(position_matrix[row][col]) == 5: #1|2|3
                    i = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5])) / 3
                elif len(position_matrix[row][col]) == 3: #1|2, 1|3, 2|3
                    i = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3])) / 2
                if i == 1.5: #1|2
                    colours.append("turquoise")
                elif i == 2:
                    colours.append("red") #1/3
                elif i == 2.5:
                    colours.append("yellow") #2/3

            elif int(position_matrix[row][col]) == 0:
                colours.append("magenta")
            elif int(position_matrix[row][col]) == 1:
                colours.append("blue")
            elif int(position_matrix[row][col]) == 2:
                colours.append("green")
            elif int(position_matrix[row][col]) == 3:
                colours.append("orange")

    center = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]


    fig = plt.figure(figsize=(15, 9))
    ax = plt.axes(projection="3d")
    ax.bar3d(x, y, z, dx, dy, dz, color=colours)
    ax.set_xlabel("True amino acid")
    ax.set_ylabel("Predicted amino acid")
    ax.set_zlabel("Number of nucleotide changes")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("FULL BAR PLOT, Z AXIS = # OF CHANGES, COLOURS = POSITION OF CHANGE")

    magenta_patch = mpatches.Patch(color='magenta', label='position of change = 0')
    blue_patch = mpatches.Patch(color='blue', label='position of change = 1')
    turq_patch = mpatches.Patch(color='turquoise', label='position of change = 1/2')
    red_patch = mpatches.Patch(color='red', label='position of change = 1/3')
    green_patch = mpatches.Patch(color='green', label='position of change = 2')
    yellow_patch = mpatches.Patch(color='yellow', label='position of change = 2/3')
    orange_patch = mpatches.Patch(color='orange', label='position of change = 3')
    plt.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch], loc='upper left', framealpha=0.5, title="Position of change")

    plt.show()

# ---------------------------------------------END---------------------------------------------------------------

# ---------------------------------------------4) HALF BAR PLOT, Z AXIS = # OF CHANGES, COLOURS = POSITION OF CHANGE---------------------------------------------------------------
if sys.argv[2] == '4':
    for row in range(20):
        for col in range(row, 20):
            x.append(row)
            y.append(col)
            z.append(0)
            dx.append(1)
            dy.append(1)
            if '|' not in codon_comparison_table[row][col]:
                dz.append(int(codon_comparison_table[row][col]))
            else:
                if len(codon_comparison_table[row][col]) == 5:
                    dz.append((int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]) + int(codon_comparison_table[row][col][4:5]))/3)
                elif len(codon_comparison_table[row][col]) == 3:
                    dz.append((int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]))/2)

            if '|' in position_matrix[row][col]:
                i = 0
                if len(position_matrix[row][col]) == 5: #1|2|3
                    i = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5])) / 3
                elif len(position_matrix[row][col]) == 3: #1|2, 1|3, 2|3
                    i = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3])) / 2
                if i == 1.5: #1|2
                    colours.append("turquoise")
                elif i == 2:
                    colours.append("red") #1/3
                elif i == 2.5:
                    colours.append("yellow") #2/3

            elif int(position_matrix[row][col]) == 0:
                colours.append("magenta")
            elif int(position_matrix[row][col]) == 1:
                colours.append("blue")
            elif int(position_matrix[row][col]) == 2:
                colours.append("green")
            elif int(position_matrix[row][col]) == 3:
                colours.append("orange")

    center = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]


    fig = plt.figure(figsize=(15, 9))
    ax = plt.axes(projection="3d")
    ax.bar3d(x, y, z, dx, dy, dz, color=colours)
    ax.set_xlabel("True amino acid")
    ax.set_ylabel("Predicted amino acid")
    ax.set_zlabel("Number of nucleotide changes")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("HALF BAR PLOT, Z AXIS = # OF CHANGES, COLOURS = POSITION OF CHANGE")

    magenta_patch = mpatches.Patch(color='magenta', label='position of change = 0')
    blue_patch = mpatches.Patch(color='blue', label='position of change = 1')
    turq_patch = mpatches.Patch(color='turquoise', label='position of change = 1/2')
    red_patch = mpatches.Patch(color='red', label='position of change = 1/3')
    green_patch = mpatches.Patch(color='green', label='position of change = 2')
    yellow_patch = mpatches.Patch(color='yellow', label='position of change = 2/3')
    orange_patch = mpatches.Patch(color='orange', label='position of change = 3')
    plt.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch], loc='upper left', framealpha=0.5, title="Position of change")

    plt.show()

# ---------------------------------------------END---------------------------------------------------------------

# ---------------------------------------------5) HALF 3D SCATTER PLOT, Z AXIS = # OF CHANGES, COLOURS = POSITION OF CHANGE---------------------------------------------------------------
if sys.argv[2] == '5':
    for row in range(20):
        for col in range(row, 20):
            x.append(row)
            y.append(col)
            # z.append(0)
            dx.append(1)
            dy.append(1)
            if '|' not in codon_comparison_table[row][col]:
                z.append(int(codon_comparison_table[row][col]))
            else:
                if len(codon_comparison_table[row][col]) == 5:
                    z.append((int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]) + int(codon_comparison_table[row][col][4:5]))/3)
                elif len(codon_comparison_table[row][col]) == 3:
                    z.append((int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]))/2)

            if '|' in position_matrix[row][col]:
                i = 0
                if len(position_matrix[row][col]) == 5: #1|2|3
                    i = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5])) / 3
                elif len(position_matrix[row][col]) == 3: #1|2, 1|3, 2|3
                    i = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3])) / 2
                if i == 1.5: #1|2
                    colours.append("turquoise")
                elif i == 2:
                    colours.append("red") #1/3
                elif i == 2.5:
                    colours.append("yellow") #2/3

            elif int(position_matrix[row][col]) == 0:
                colours.append("magenta")
            elif int(position_matrix[row][col]) == 1:
                colours.append("blue")
            elif int(position_matrix[row][col]) == 2:
                colours.append("green")
            elif int(position_matrix[row][col]) == 3:
                colours.append("orange")

    center = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]


    fig = plt.figure(figsize=(15, 9))
    ax = plt.axes(projection="3d")
    ax.scatter(x, y, z, c=colours)
    ax.set_xlabel("True amino acid")
    ax.set_ylabel("Predicted amino acid")
    ax.set_zlabel("Number of nucleotide changes")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("HALF 3D SCATTER PLOT, Z AXIS = # OF CHANGES, COLOURS = POSITION OF CHANGE")

    magenta_patch = mpatches.Patch(color='magenta', label='position of change = 0')
    blue_patch = mpatches.Patch(color='blue', label='position of change = 1')
    turq_patch = mpatches.Patch(color='turquoise', label='position of change = 1/2')
    red_patch = mpatches.Patch(color='red', label='position of change = 1/3')
    green_patch = mpatches.Patch(color='green', label='position of change = 2')
    yellow_patch = mpatches.Patch(color='yellow', label='position of change = 2/3')
    orange_patch = mpatches.Patch(color='orange', label='position of change = 3')
    plt.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch], loc='upper left', framealpha=0.5, title="Position of change")

    plt.show()

# ---------------------------------------------END---------------------------------------------------------------

# ---------------------------------------------6) HALF 3D SCATTER PLOT, Z AXIS = POSITION OF CHANGE, COLOURS = # OF CHANGES---------------------------------------------------------------
if sys.argv[2] == '6':
    for row in range(20):
        for col in range(row, 20):
            x.append(row)
            y.append(col)
            # z.append(0)
            dx.append(1)
            dy.append(1)
            if '|' not in position_matrix[row][col]:
                z.append(int(position_matrix[row][col]))
            else:
                if len(position_matrix[row][col]) == 5:
                    z.append((int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5]))/3)
                elif len(position_matrix[row][col]) == 3:
                    z.append((int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]))/2)

            if '|' in codon_comparison_table[row][col]:
                i = 0
                if len(codon_comparison_table[row][col]) == 5: #1|2|3
                    i = (int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]) + int(codon_comparison_table[row][col][4:5])) / 3
                elif len(codon_comparison_table[row][col]) == 3: #1|2, 1|3, 2|3
                    i = (int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3])) / 2
                if i == 1.5: #1|2
                    colours.append("turquoise")
                elif i == 2:
                    colours.append("red") #1/3
                elif i == 2.5:
                    colours.append("yellow") #2/3

            elif int(codon_comparison_table[row][col]) == 0:
                colours.append("magenta")
            elif int(codon_comparison_table[row][col]) == 1:
                colours.append("blue")
            elif int(codon_comparison_table[row][col]) == 2:
                colours.append("green")
            elif int(codon_comparison_table[row][col]) == 3:
                colours.append("orange")

    center = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]


    fig = plt.figure(figsize=(15, 9))
    ax = plt.axes(projection="3d")
    ax.scatter(x, y, z, c=colours)
    ax.set_xlabel("True amino acid")
    ax.set_ylabel("Predicted amino acid")
    ax.set_zlabel("Position of change")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("HALF 3D SCATTER PLOT, Z AXIS = POSITION OF CHANGE, COLOURS = # OF CHANGES")

    magenta_patch = mpatches.Patch(color='magenta', label='# of changes = 0')
    blue_patch = mpatches.Patch(color='blue', label='# of changes = 1')
    turq_patch = mpatches.Patch(color='turquoise', label='# of changes = 1/2')
    red_patch = mpatches.Patch(color='red', label='# of changes = 1/3')
    green_patch = mpatches.Patch(color='green', label='# of changes = 2')
    yellow_patch = mpatches.Patch(color='yellow', label='# of changes = 2/3')
    orange_patch = mpatches.Patch(color='orange', label='# of changes = 3')
    plt.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch], loc='upper left', framealpha=0.5, title="Number of changes")

    plt.show()

# ---------------------------------------------END---------------------------------------------------------------

# ---------------------------------------------7) FULL 2D SCATTER PLOT, COLOURS = POSITION OF CHANGE, SIZES = # OF CHANGES---------------------------------------------------------------
if sys.argv[2] == '7':
    c = []
    s = []
    edges = []
    pos = 0
    highest_freq = {0:18, 1:2, 2:4, 3:4, 4:2, 5:4, 6:5, 7:8, 8:4, 9:4, 10:4, 11:4, 12:13, 13:14, 14:13, 15:16, 16:14, 17:18, 18:17, 19:17} #row:col


    for row in range(20):
        for col in range(20):
            x.append(row)
            y.append(col)
            dx.append(1)
            dy.append(1)
            if '|' not in position_matrix[row][col]:
                pos = int(position_matrix[row][col])
                if pos == 0:
                    c.append("grey")
                    edges.append("grey")
                elif pos == 1:
                    c.append("royalblue")
                    edges.append("royalblue")
                elif pos == 2:
                    c.append("green")
                    edges.append("green")
                elif pos == 3:
                    c.append("orange")
                    edges.append("orange")
            else:
                if len(position_matrix[row][col]) == 5:
                    # pos = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]) + int(position_matrix[row][col][4:5]))/3
                    c.append("hotpink")
                    edges.append("hotpink")
                elif len(position_matrix[row][col]) == 3:
                    pos = (int(position_matrix[row][col][0:1]) + int(position_matrix[row][col][2:3]))/2
                    if pos == 1.5:
                        c.append("turquoise")
                        edges.append("turquoise")
                    elif pos == 2:
                        c.append("red")
                        edges.append("red")
                    elif pos == 2.5:
                        c.append("yellow")
                        edges.append("yellow")
                    
            if '|' in codon_comparison_table[row][col]:
                i = 0
                if len(codon_comparison_table[row][col]) == 5: #1|2|3
                    i = (int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3]) + int(codon_comparison_table[row][col][4:5])) / 3
                elif len(codon_comparison_table[row][col]) == 3: #1|2, 1|3, 2|3
                    i = (int(codon_comparison_table[row][col][0:1]) + int(codon_comparison_table[row][col][2:3])) / 2
                if i == 1.5: #1|2
                    s.append(1.5*120)
                elif i == 2:
                    s.append(2.01*120) #1/3
                elif i == 2.5:
                    s.append(2.5*120) #2/3

            elif int(codon_comparison_table[row][col]) == 0:
                s.append(1*40)
            elif int(codon_comparison_table[row][col]) == 1:
                s.append(1*120)
            elif int(codon_comparison_table[row][col]) == 2:
                s.append(2*120)
            elif int(codon_comparison_table[row][col]) == 3:
                s.append(3*130)
            
            if (row, col) in highest_freq.items():
                pos = len(edges)
                edges[pos - 1] = "black"

    center = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

    fig, ax = plt.subplots(figsize=(15, 9))
    scatter = ax.scatter(y, x, c=c, s=s, edgecolors=edges, linewidth=3)

    handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    labels = ["None", '1', '1 or 2', '1 or 3', '2', '2 or 3', '3']
    legend2 = ax.legend(handles, labels, bbox_to_anchor=(1.1, 0.7), loc="upper right", framealpha=0.5, title="# of changes")
    ax.add_artist(legend2)

    ax.set_xlabel("Predicted amino acid")
    ax.set_ylabel("True amino acid")
    ax.set_xticks(center, amino_acids)
    ax.set_yticks(center, amino_acids)
    plt.title("FULL 2D SCATTER PLOT, COLOURS = POSITION OF CHANGE, SIZES = # OF CHANGES")

    magenta_patch = mpatches.Patch(color='grey', label='No change')
    blue_patch = mpatches.Patch(color='royalblue', label='1')
    turq_patch = mpatches.Patch(color='turquoise', label='1/2')
    red_patch = mpatches.Patch(color='red', label='1/3')
    green_patch = mpatches.Patch(color='green', label='2')
    yellow_patch = mpatches.Patch(color='yellow', label='2/3')
    orange_patch = mpatches.Patch(color='orange', label='3')
    hotpink_patch = mpatches.Patch(color='hotpink', label='1/2/3')
    legend1 = ax.legend(handles=[magenta_patch, blue_patch, turq_patch, red_patch, green_patch, yellow_patch, orange_patch, hotpink_patch], bbox_to_anchor=(1.13, 1), loc='upper right', framealpha=0.5, title="Position of change")
    ax.add_artist(legend1)

    plt.show()

    print("No change: " + str(c.count('grey')))
    print("Position 1: " + str(c.count('royalblue')))
    print("Position 1 or 2: " + str(c.count('turquoise')))
    print("Position 1 or 3: " + str(c.count('red')))
    print("Position 2: " + str(c.count('green')))
    print("Position 2 or 3: " + str(c.count('yellow')))
    print("Position 3: " + str(c.count('orange')))
    print("Position 1, 2 or 3: " + str(c.count('hotpink')))

# ---------------------------------------------END---------------------------------------------------------------