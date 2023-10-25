from cmath import sqrt
from operator import index
import random
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import matplotlib.cm as cm
import numpy
from scipy.stats import kstest
from scipy.stats import normaltest
from scipy.stats import shapiro 
import sklearn.metrics as metrics

#generate normalised substitution matrix list from .csv
normalised_substitution_matrix = list(csv.reader(open("/home/larajuneb/Honours/PROJECT_(721)/Coding/BioinfoHonours/oversampled_normalised_matrix.csv")))
normalised_substitution_matrix.remove(normalised_substitution_matrix[0])
for row in range(len(normalised_substitution_matrix)):
    normalised_substitution_matrix[row].remove(normalised_substitution_matrix[row][0])

#SeqPredNN codon matrix
SeqPredNN_codon_matrix = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix = [[0 for x in range(61)] for y in range(61)]

SeqPredNN_codon_matrix_PLOTS = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix_PLOTS = [[0 for x in range(61)] for y in range(61)]
# Higgs_codon_matrix = [[0 for x in range(61)] for y in range(61)]

SeqPredNN_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]

nucleotides = ['U', 'C', 'A', 'G']
codons = []
codons_excl_stop = []
codons_leftover = []
codon = ""
stop = ['UAA', 'UAG', 'UGA']

#generate all codons
for first in nucleotides:
    for second in nucleotides:
        for third in nucleotides:
            codon = first + second + third
            codons.append(codon)
            codons_leftover.append(codon)
            if codon not in stop:
                codons_excl_stop.append(codon)

number_of_codons_per_aa = {"HIS": 2, "ARG": 6, "LYS": 2, "GLN": 2, "GLU": 2, "ASP": 2, "ASN": 2, "GLY": 4, "ALA": 4, "SER": 6, "THR": 4, "PRO": 4, "CYS": 2, "VAL": 4, "ILE": 3, "MET": 1, "LEU": 6, "PHE": 2, "TYR": 2, "TRP": 1, "stop": 3} #MAY NOT BE NECESSARY, CAN JUST GET LENGTH OF AMINO ACID CODON ARRAY BELOW

amino_acids = ["HIS", "ARG", "LYS", "GLN", "GLU", "ASP", "ASN", "GLY", "ALA", "SER", "THR", "PRO", "CYS", "VAL", "ILE", "MET", "LEU", "PHE", "TYR", "TRP"]

codons_per_aa = {"HIS": ['CAU', 'CAC'], "ARG": ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], "LYS": ['AAA', 'AAG'], "GLN": ['CAA', 'CAG'], "GLU": ['GAA', 'GAG'], "ASP": ['GAU', 'GAC'], "ASN": ['AAU', 'AAC'], "GLY": ['GGU', 'GGC', 'GGA', 'GGG'], "ALA": ['GCU', 'GCC', 'GCA', 'GCG'], "SER": ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], "THR": ['ACU', 'ACC', 'ACA', 'ACG'], "PRO": ['CCU', 'CCC', 'CCA', 'CCG'], "CYS": ['UGU', 'UGC'], "VAL": ['GUU', 'GUC', 'GUA', 'GUG'], "ILE": ['AUU', 'AUC', 'AUA'], "MET": ['AUG'], "LEU": ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'], "PHE": ['UUU', 'UUC'], "TYR": ['UAU', 'UAC'], "TRP": ['UGG'], "stop": ['UAA', 'UAG', 'UGA']}


N = 5.7 #using the Freeland Hurst weights below
FreelandHurst_mutation_weights = {}
FreelandHurst_mutation_weights["3"] = 1/N
FreelandHurst_mutation_weights["1ts"] = 1/N
FreelandHurst_mutation_weights["1tv"] = 0.5/N
FreelandHurst_mutation_weights["2ts"] = 0.5/N
FreelandHurst_mutation_weights["2tv"] = 0.1/N
FreelandHurst_mutation_weights["otherwise"] = 0 #more than 1 position that differs

polar_requirement_scale = {"HIS": 8.4, "ARG": 9.1, "LYS": 10.1, "GLN": 8.6, "GLU": 12.5, "ASP": 13, "ASN": 10, "GLY": 7.9, "ALA": 7, "SER": 7.5, "THR": 6.6, "PRO": 6.6, "CYS": 11.5, "VAL": 5.6, "ILE": 4.9, "MET": 5.3, "LEU": 4.9, "PHE": 5, "TYR": 5.7, "TRP": 5.3}

all_N = []
all_N_codons = []

purines = ['A', 'G']
pyrimidines = ['C', 'U']

#identify if codon mutation is a transition or transversion
def tv_or_ts(true, mutant):
    if (true in purines and mutant in purines) or (true in pyrimidines and mutant in pyrimidines):
        return "ts"
    elif (true in purines and mutant in pyrimidines) or (true in pyrimidines and mutant in purines):
        return "tv"

#calculate N constant for Freeland and Hurst weights
def calculate_N():
    for true_codon in codons:
        sum_for_N = 0
        first = 0
        second = 0
        third = 0
        # print
        for mutant_codon in codons:
            if true_codon != mutant_codon:
                sum_of_diffs = 0
                first = 0
                second = 0
                third = 0
                if true_codon[0:1] != mutant_codon[0:1]:
                    first = 1
                if true_codon[1:2] != mutant_codon[1:2]:
                    second = 1
                if true_codon[2:3] != mutant_codon[2:3]:
                    third = 1
                sum_of_diffs = first + second + third

                if sum_of_diffs <= 1: #weight for more than one position difference is 0, so don't need to add that to sum_for_N
                    if first == 1: #differ in first position only
                        if tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "ts":
                            sum_for_N += FreelandHurst_mutation_weights.get("1ts")
                        elif tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "tv":
                            sum_for_N += FreelandHurst_mutation_weights.get("1tv")
                    elif second == 1: #differ in second position only
                        if tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "ts":
                            sum_for_N += FreelandHurst_mutation_weights.get("2ts")
                        elif tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "tv":
                            sum_for_N += FreelandHurst_mutation_weights.get("2tv")
                    elif third == 1: #differ in third position only
                        sum_for_N += FreelandHurst_mutation_weights.get("3")
        all_N.append(sum_for_N)
        all_N_codons.append(true_codon)

#calculate probability of codon mutation using Freeland and Hurst mutation weights
def get_codon_mutation_prob(true_codon, mutant_codon):
    if true_codon != mutant_codon:
        sum_of_diffs = 0
        first = 0
        second = 0
        third = 0
        if true_codon[0:1] != mutant_codon[0:1]:
            first = 1
        if true_codon[1:2] != mutant_codon[1:2]:
            second = 1
        if true_codon[2:3] != mutant_codon[2:3]:
            third = 1
        sum_of_diffs = first + second + third

        if sum_of_diffs <= 1:# differ at only 1 position
            if first == 1: #differ in first position only
                if tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "ts":
                    return(FreelandHurst_mutation_weights.get("1ts"))
                elif tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "tv":
                    return(FreelandHurst_mutation_weights.get("1tv"))
            elif second == 1: #differ in second position only
                if tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "ts":
                    return(FreelandHurst_mutation_weights.get("2ts"))
                elif tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "tv":
                    return(FreelandHurst_mutation_weights.get("2tv"))
            elif third == 1: #differ in third position only
                return(FreelandHurst_mutation_weights.get("3"))
        elif sum_of_diffs > 1: #differ at more than one position
            return '-'
    elif true_codon == mutant_codon:
        return '*'

#get the amino acid that a codon codes for
def get_aa_for_codon(codon, codon_dict):
    for key, value in codon_dict.items():
        if codon in value:
            return key

#calculate the cost of a codon mutation
def get_cost(true_codon_index, mutant_codon_index, mode, codon_dict, codons_without_stop, plot):
    #get codon string name from it's index
    true_codon = codons_without_stop[true_codon_index]
    mutant_codon = codons_without_stop[mutant_codon_index]
    #get amino acid name from codon name
    true_aa = get_aa_for_codon(true_codon, codon_dict)
    mutant_aa = get_aa_for_codon(mutant_codon, codon_dict)
    #get amino acid name from amino acid index
    true_aa_index = amino_acids.index(true_aa)
    mutant_aa_index = amino_acids.index(mutant_aa)

    aa_difference = 0
    if mode == "SeqPredNN":
        aa_difference = 1 - float(normalised_substitution_matrix[true_aa_index][mutant_aa_index])
    elif mode == "Koonin":
        aa_difference = pow((polar_requirement_scale.get(true_aa) - polar_requirement_scale.get(mutant_aa)), 2)

    codon_mutation_prob = get_codon_mutation_prob(true_codon, mutant_codon)
    if codon_mutation_prob == '*':
        cost = '*'
    elif codon_mutation_prob == '-':
        cost = '-'
    else:
        cost = float(codon_mutation_prob) * float(aa_difference)
    if plot == True and isinstance(cost, str):
        cost = -10
    return (cost)

#normalise the cost value to be between 0 and 100
def normalise_cost(min, max, value):
    z = ((value - min) / (max - min)) * 100
    return(z)

#calculate the overall cost of the code
def get_code_cost(cost_array):
    code_cost = 0
    for row in range(61):
        for cell in range(61):
            if not isinstance(cost_array[row][cell], str):
                code_cost += cost_array[row][cell]
    return code_cost

def regenerate_codon_matrix(array):
    for first in nucleotides:
        for second in nucleotides:
            for third in nucleotides:
                codon = first + second + third
                array.append(codon)
    return(array)

#check if a list of costs in normally distributed using the Kolmogorov-Smirnov Test
def is_it_gaussian(dataset):
    kstest(dataset, 'norm')
    return 0
#generate 10,000 random assignments of codons and calculate costs for each
def generate_sample_set(sample_size, SeqPredNN_sample_code_costs, SeqPredNN_sample_code_costs_NORM, Koonin_sample_code_costs, Koonin_sample_code_costs_NORM):
    random_codon_assignments = {} #similar to codons_per_aa dict, but instead of true codon assignments, the codons are assigned randomly to amino acids
    leftover = []
    stop = []
    code_costs = []
    norm_code_costs = []
    no_stop_codons = []
    for j in range(2):
        for i in range(sample_size):
            temp_cost_matrix = [[0 for x in range(61)] for y in range(61)]
            norm_cost_matrix = [[0 for x in range(61)] for y in range(61)]
            matrix_min_max_check = []
            minimum = 0
            maximum = 0
            leftover = regenerate_codon_matrix(leftover)
            no_stop_codons = regenerate_codon_matrix(no_stop_codons)
            for key, value in number_of_codons_per_aa.items():
                random_codon_assignments[key] = [] #a new random codon assignment for each sample, key = amino acid, value = array of codons
                for i in range(number_of_codons_per_aa[key]): #loop through as many times as there are codons for that amino acid
                    codon = random.choice(leftover) #choose a random codon for the amino acid 
                    leftover.remove(codon) #remove that codon from the list of available codons
                    random_codon_assignments[key].append(codon) #add the randomly chose codon to the list of codons for that amino acid
            stop = random_codon_assignments.get("stop")
            #random codon assignment have been made, now calculate cost matrix and overall cost
            for codon in stop:
                no_stop_codons.remove(codon)
            
            for row in range(61):
                for cell in range(61):
                    if j == 0:
                        temp_cost_matrix[row][cell] = get_cost(row, cell, "SeqPredNN", random_codon_assignments, no_stop_codons, False)
                    elif j == 1:
                        temp_cost_matrix[row][cell] = get_cost(row, cell, "Koonin", random_codon_assignments, no_stop_codons, False)

                    if not isinstance(temp_cost_matrix[row][cell], str):
                        matrix_min_max_check.append(temp_cost_matrix[row][cell])
            #get cosde cost for raw data matrix
            minimum = min(matrix_min_max_check)
            maximum = max(matrix_min_max_check)
            code_costs.append(get_code_cost(temp_cost_matrix))

            #get code cost for normalised data matrix
            norm_cost_matrix = normalise_matrix(minimum, maximum, temp_cost_matrix)

            norm_code_costs.append(get_code_cost(norm_cost_matrix))

            if j == 0:
                SeqPredNN_sample_code_costs.append(get_code_cost(temp_cost_matrix))
                SeqPredNN_sample_code_costs_NORM.append(get_code_cost(norm_cost_matrix))
            if j == 1:
                Koonin_sample_code_costs.append(get_code_cost(temp_cost_matrix))
                Koonin_sample_code_costs_NORM.append(get_code_cost(norm_cost_matrix))

            random_codon_assignments.clear() #clear random codon assignment for next sample
            no_stop_codons.clear()
            stop.clear()
            temp_cost_matrix = [[0 for x in range(61)] for y in range(61)]
            norm_cost_matrix = [[0 for x in range(61)] for y in range(61)]

        if j == 0:
            print("************************************ SEQPREDNN ************************************")
        elif j == 1:
            print("************************************ KOONIN ************************************")
        print("----------------------- series -----------------------")
        series = pd.Series(code_costs)
        print(series.describe())
        print("----------------------- norm series -----------------------")
        series_norm = pd.Series(norm_code_costs)
        print(series_norm.describe())
        
        if j == 0:
            plot_samples(sample_size, code_costs, "SeqPredNN", False)
            plot_samples(sample_size, norm_code_costs, "SeqPredNN", True)
        if j == 1:
            plot_samples(sample_size, code_costs, "Koonin", False)
            plot_samples(sample_size, norm_code_costs, "Koonin", True)

        code_costs.clear()
        norm_code_costs.clear()
    
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("SeqPredNN_sample_code_costs: " + str(shapiro(SeqPredNN_sample_code_costs)))
    if shapiro(SeqPredNN_sample_code_costs).pvalue > 0.05:
        print("Normal = TRUE")
    print("SeqPredNN_sample_code_costs_NORM: " + str(shapiro(SeqPredNN_sample_code_costs_NORM)))
    if shapiro(SeqPredNN_sample_code_costs_NORM).pvalue > 0.05:
        print("Normal = TRUE")
    print("Koonin_sample_code_costs: " + str(shapiro(Koonin_sample_code_costs)))
    if shapiro(Koonin_sample_code_costs).pvalue > 0.05:
        print("Normal = TRUE")
    print("Koonin_sample_code_costs_NORM: " + str(shapiro(Koonin_sample_code_costs_NORM)))
    if shapiro(Koonin_sample_code_costs_NORM).pvalue > 0.05:
        print("Normal = TRUE")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("SeqPredNN_sample_code_costs: " + str(normaltest(SeqPredNN_sample_code_costs)))
    if normaltest(SeqPredNN_sample_code_costs).pvalue > 0.05:
        print("Normal = TRUE")
    print("SeqPredNN_sample_code_costs_NORM: " + str(normaltest(SeqPredNN_sample_code_costs_NORM)))
    if normaltest(SeqPredNN_sample_code_costs_NORM).pvalue > 0.05:
        print("Normal = TRUE")
    print("Koonin_sample_code_costs: " + str(normaltest(Koonin_sample_code_costs)))
    if normaltest(Koonin_sample_code_costs).pvalue > 0.05:
        print("Normal = TRUE")
    print("Koonin_sample_code_costs_NORM: " + str(normaltest(Koonin_sample_code_costs_NORM)))
    if normaltest(Koonin_sample_code_costs_NORM).pvalue > 0.05:
        print("Normal = TRUE")
    
#plot bar graphs for samples produced
def plot_samples(sample_size, costs, model, normalised):
    costs_temp = []
    for item in costs:
        costs_temp.append(item)
    num_bins = math.isqrt(sample_size)
    mini = float(min(costs))
    maxi = float(max(costs))
    step = (maxi - mini) / num_bins
    bins = {}
    bin_counts = {}
    val = mini
    x = []
    x_bins = []
    y = []
    for i in range(num_bins):
        if i == (num_bins - 1):
            bins[i] = [val, maxi]
        else:
            bins[i] = [val, val + step]
        bin_counts[i] = 0
        val += step

    count = 0
    for key, value in bins.items():
        for item in costs_temp:
            if item >= value[0] and item <= value[1]:
                bin_counts[key] += 1
        x_bins.append("[" + str(round(bins[key][0], 3)) + " - " + str(round(bins[key][1], 3)) + "]")
        x.append(count)
        y.append(bin_counts[key])
        count += 1

    plt.figure(figsize=(10, 8))
    plt.bar(x_bins, y)
    plt.xticks(rotation = 20)
    if normalised == False:
        plt.title("Number of occurrences of randomly simulated code costs per bracket for models generated\n using the " + model + " method")
    if normalised == True:
        plt.title("Number of occurrences of randomly simulated code costs per bracket for NORMALISED models\n generated using the " + model + " method")
    plt.xlabel("Code cost")
    plt.ylabel("Number of occurrences")
    plt.show()

    return 0

#make csv files of all matrices
def make_csvs():
    filename = "SeqPredNN_cost_matrix.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in SeqPredNN_codon_matrix:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "Koonin_cost_matrix.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in Koonin_codon_matrix:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "SeqPredNN_cost_matrix_NORM.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in SeqPredNN_codon_matrix_NORM:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "Koonin_cost_matrix_NORM.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in Koonin_codon_matrix_NORM:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

#start cost calculations
def calculate_SeqPredNN_and_Koonin_matrices():
    for row in range(61):
        for cell in range(61):
            SeqPredNN_codon_matrix[row][cell] = get_cost(row, cell, "SeqPredNN", codons_per_aa, codons_excl_stop, False)
            Koonin_codon_matrix[row][cell] = get_cost(row, cell, "Koonin", codons_per_aa, codons_excl_stop, False)

            if not isinstance(SeqPredNN_codon_matrix[row][cell], str):
                SeqPredNN_check.append(SeqPredNN_codon_matrix[row][cell])

            if not isinstance(Koonin_codon_matrix[row][cell], str):
                Koonin_check.append(Koonin_codon_matrix[row][cell])

            SeqPredNN_codon_matrix_PLOTS[row][cell] = get_cost(row, cell, "SeqPredNN", codons_per_aa, codons_excl_stop, True)
            Koonin_codon_matrix_PLOTS[row][cell] = get_cost(row, cell, "Koonin", codons_per_aa, codons_excl_stop, True)

def normalise_matrix(minimum, maximum, original_matrix):
    normalised = [[0 for x in range(61)] for y in range(61)]

    for row in range(len(original_matrix)):
        for cell in range(len(original_matrix)):
            if not isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
            if isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = original_matrix[row][cell]

    return(normalised)

# def confusion_matrix(self, true_residues, predicted_residues, normalize, file_name):
def confusion_matrix(plot_matrix, original_matrix, minimum, maximum, filename):

    fig, ax = plt.subplots(figsize=(42, 35))
    plt.rcParams['font.size'] = 11
    img = ax.imshow(plot_matrix)
    ax.set_xticks(np.arange(61))
    ax.set_yticks(np.arange(61))
    ax.set_xticklabels(codons_excl_stop, rotation = 90)
    ax.set_yticklabels(codons_excl_stop)
    ax.tick_params(labelsize=30, pad=14, length=14, width=3)
    ax.set_xlabel('Mutant codon', fontsize=70, labelpad=24)
    ax.set_ylabel('True codon', fontsize=70, labelpad=24)

    # psm = ax.pcolormesh(plot_matrix, cmap='viridis', vmin=minimum, vmax=maximum)
    # fig.colorbar(psm, ax=ax)
    for i in range(61):
        for j in range(61):
            if not isinstance(original_matrix[i][j], str):
                ax.text(j, i, round(original_matrix[i][j], 3), ha="center", va="center", color="black")
            else:
                ax.text(j, i, original_matrix[i][j], ha="center", va="center", color="black")

    cbar = fig.colorbar(img, ax=ax, pad=0.01)
    cbar.ax.tick_params(labelsize=63, pad=14, length=14, width=3)
    fig.tight_layout()
    plt.savefig(f'{filename}.png')
    plt.close()

    
    # plt.close()

SeqPredNN_cost_min = 0
SeqPredNN_cost_max = 0
Koonin_cost_min = 0
Koonin_cost_max = 0

SeqPredNN_check = []
Koonin_check = []

#run calculations for SeqPredNN and Koonin matrices
calculate_SeqPredNN_and_Koonin_matrices()

#calculate min and max for both matrices
SeqPredNN_cost_min = min(SeqPredNN_check)
SeqPredNN_cost_max = max(SeqPredNN_check)
Koonin_cost_min = min(Koonin_check)
Koonin_cost_max = max(Koonin_check)

#create normalised codon mutation cost matrices
SeqPredNN_codon_matrix_NORM = normalise_matrix(SeqPredNN_cost_min, SeqPredNN_cost_max, SeqPredNN_codon_matrix)
Koonin_codon_matrix_NORM = normalise_matrix(Koonin_cost_min, Koonin_cost_max, Koonin_codon_matrix)

#get overall code cost for SeqPredNN and Koonin matrices
SeqPredNN_code_cost = get_code_cost(SeqPredNN_codon_matrix)
Koonin_code_cost = get_code_cost(Koonin_codon_matrix)

print("SeqPredNN code cost: " + str(SeqPredNN_code_cost))
print("SeqPredNN min: " + str(min(SeqPredNN_check)))
print("SeqPredNN max: " + str(max(SeqPredNN_check)))
print("Koonin code cost: " + str(Koonin_code_cost))
print("Koonin min: " + str(min(Koonin_check)))
print("Koonin max: " + str(max(Koonin_check)))

print("NORM SeqPredNN code cost: " + str(get_code_cost(SeqPredNN_codon_matrix_NORM)))
print("NORM Koonin code cost: " + str(get_code_cost(Koonin_codon_matrix_NORM)))

#generate 10,000 random codon assignments for SeqPredNN model and calculate code costs for each assignment
SeqPredNN_sample_code_costs = []
SeqPredNN_sample_code_costs_NORM = []
Koonin_sample_code_costs = []
Koonin_sample_code_costs_NORM = []
# generate_sample_set(800, SeqPredNN_sample_code_costs, SeqPredNN_sample_code_costs_NORM, Koonin_sample_code_costs, Koonin_sample_code_costs_NORM)

#generate .csv files
make_csvs()

#generate plots
confusion_matrix(SeqPredNN_codon_matrix_PLOTS, SeqPredNN_codon_matrix, min(SeqPredNN_check), max(SeqPredNN_check), "SeqPredNN_codon_matrix")
confusion_matrix(Koonin_codon_matrix_PLOTS, Koonin_codon_matrix, min(Koonin_check), max(Koonin_check), "Koonin_codon_matrix")

symmetrical_S = []
symmetrical_K = []
for i in range(61):
    for j in range(61):
        if SeqPredNN_codon_matrix[i][j] == SeqPredNN_codon_matrix[j][i]:
            symmetrical_S.append(True)
        elif SeqPredNN_codon_matrix[i][j] != SeqPredNN_codon_matrix[j][i]:
            symmetrical_S.append(False)
            print("SEQPREDNN NOT SYMMETRICAL: i=" + str(i) + " j=" + str(j))
            
        if Koonin_codon_matrix[i][j] == Koonin_codon_matrix[j][i]:
            symmetrical_K.append(True)
        elif Koonin_codon_matrix[i][j] != Koonin_codon_matrix[j][i]:
            symmetrical_K.append(False)
            print("KOONIN NOT SYMMETRICAL: i=" + str(i) + " j=" + str(j))

if False in symmetrical_S:
    print("SEQPREDNN NOT SYMMETRICAL")

if False in symmetrical_K:
    print("KOONIN NOT SYMMETRICAL")
