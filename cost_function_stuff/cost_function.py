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
from scipy.stats import spearmanr 
from scipy.stats import ks_2samp
import sklearn.metrics as metrics

#generate normalised substitution matrix list from .csv
normalised_substitution_matrix = list(csv.reader(open("/home/larajuneb/Honours/PROJECT_(721)/Coding/BioinfoHonours/oversampled_normalised_matrix.csv")))
normalised_substitution_matrix.remove(normalised_substitution_matrix[0])


Higgs_distance_matrix = list(csv.reader(open("Higgs_aa_distance_matrix.csv")))
Higgs_distance_matrix.remove(Higgs_distance_matrix[0])

for row in range(len(normalised_substitution_matrix)):
    normalised_substitution_matrix[row].remove(normalised_substitution_matrix[row][0])

for row in range(len(Higgs_distance_matrix)):
    Higgs_distance_matrix[row].remove(Higgs_distance_matrix[row][0])

#SeqPredNN codon matrix
SeqPredNN_codon_matrix = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix = [[0 for x in range(61)] for y in range(61)]
Higgs_codon_matrix = [[0 for x in range(61)] for y in range(61)]

neutral_substitution_matrix = [[0 for x in range(61)] for y in range(61)]
neutral_substitution_matrix_PLOTS = [[0 for x in range(61)] for y in range(61)]

SeqPredNN_codon_matrix_PLOTS = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix_PLOTS = [[0 for x in range(61)] for y in range(61)]
Higgs_codon_matrix_PLOTS = [[0 for x in range(61)] for y in range(61)]

SeqPredNN_codon_matrix_NORM_PLOTS = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix_NORM_PLOTS = [[0 for x in range(61)] for y in range(61)]
Higgs_codon_matrix_NORM_PLOTS = [[0 for x in range(61)] for y in range(61)]

SeqPredNN_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
Higgs_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]


SeqPredNN_primordial_codon_matrix = [[0 for x in range(16)] for y in range(16)]
Koonin_primordial_codon_matrix = [[0 for x in range(16)] for y in range(16)]
Higgs_primordial_codon_matrix = [[0 for x in range(16)] for y in range(16)]
SeqPredNN_primordial_codon_matrix_NORM = [[0 for x in range(16)] for y in range(16)]
Koonin_primordial_codon_matrix_NORM = [[0 for x in range(16)] for y in range(16)]
Higgs_primordial_codon_matrix_NORM = [[0 for x in range(16)] for y in range(16)]
SeqPredNN_primordial_check = []
Koonin_primordial_check = []
Higgs_primordial_check = []

SeqPredNN_check = []
SeqPredNN_check_NORM = []
Koonin_check = []
Koonin_check_NORM = []
neutral_check = []
Higgs_check = []
Higgs_check_NORM = []

SeqPredNN_sample_code_costs = []
SeqPredNN_sample_code_costs_NORM = []
Koonin_sample_code_costs = []
Koonin_sample_code_costs_NORM = []
Higgs_sample_code_costs = []
Higgs_sample_code_costs_NORM = []

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

Higgs_amino_acid_order = ['PHE','LEU','ILE','MET','VAL','SER','PRO','THR','ALA','TYR','HIS','GLN','ASN','LYS','ASP','GLU','CYS','TRP','ARG','GLY']

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
    if mode != "Higgs":
        true_aa_index = amino_acids.index(true_aa)
        mutant_aa_index = amino_acids.index(mutant_aa)
    elif mode == "Higgs":
        true_aa_index = Higgs_amino_acid_order.index(true_aa)
        mutant_aa_index = Higgs_amino_acid_order.index(mutant_aa)

    aa_difference = 0
    if mode == "SeqPredNN":
        aa_difference = 1 - float(normalised_substitution_matrix[true_aa_index][mutant_aa_index])
    elif mode == "Koonin":
        aa_difference = pow((polar_requirement_scale.get(true_aa) - polar_requirement_scale.get(mutant_aa)), 2)
    elif mode == "Higgs":
        aa_difference = Higgs_distance_matrix[true_aa_index][mutant_aa_index]
    elif mode == "neutral":
        aa_difference = 1

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

def generate_sample_set(sample_size, sample_code_costs, sample_code_costs_NORM, mode_code_cost, mode_code_cost_NORM, mode):
    random_codon_assignments = {} #similar to codons_per_aa dict, but instead of true codon assignments, the codons are assigned randomly to amino acids
    leftover = []
    stop = []
    code_costs = []
    norm_code_costs = []
    no_stop_codons = []
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
                temp_cost_matrix[row][cell] = get_cost(row, cell, mode, random_codon_assignments, no_stop_codons, False)

                if not isinstance(temp_cost_matrix[row][cell], str):
                    matrix_min_max_check.append(temp_cost_matrix[row][cell])
        #get code cost for raw data matrix
        minimum = min(matrix_min_max_check)
        maximum = max(matrix_min_max_check)
        code_costs.append(get_code_cost(temp_cost_matrix))

        #get code cost for normalised data matrix
        norm_cost_matrix = normalise_matrix_SAMPLES(minimum, maximum, temp_cost_matrix)

        norm_code_costs.append(get_code_cost(norm_cost_matrix))

        sample_code_costs.append(get_code_cost(temp_cost_matrix))
        sample_code_costs_NORM.append(get_code_cost(norm_cost_matrix))

        random_codon_assignments.clear() #clear random codon assignment for next sample
        no_stop_codons.clear()
        stop.clear()
        temp_cost_matrix = [[0 for x in range(61)] for y in range(61)]
        norm_cost_matrix = [[0 for x in range(61)] for y in range(61)]
        
    sample_set_stats(code_costs, norm_code_costs, mode)

    plot_samples(sample_size, code_costs, mode, False, mode_code_cost)
    plot_samples(sample_size, norm_code_costs, mode, True, mode_code_cost_NORM)

    code_costs.clear()
    norm_code_costs.clear()

def sample_set_stats(costs, norm_costs, mode):
    series = pd.Series(costs)
    series_norm = pd.Series(norm_costs)

    filename = "stats/" + mode + "_code_cost_samples_stats.csv"
    with open(filename, mode="w") as file:
        file.write("mode, count, mean, std, min, 25%, 50%, 75%, max\n")

        file.write(mode + "," + str(series.describe()[0]) + "," + str(series.describe()[1]) + "," + str(series.describe()[2]) + "," + str(series.describe()[3]) + "," + str(series.describe()[4]) + "," + str(series.describe()[5]) + "," + str(series.describe()[6]) + "," + str(series.describe()[7]) + "\n")

        file.write(mode + "NORM," + str(series_norm.describe()[0]) + "," + str(series_norm.describe()[1]) + "," + str(series_norm.describe()[2]) + "," + str(series_norm.describe()[3]) + "," + str(series_norm.describe()[4]) + "," + str(series_norm.describe()[5]) + "," + str(series_norm.describe()[6]) + "," + str(series_norm.describe()[7]) + "\n")

#plot bar graphs for samples produced
def plot_samples(sample_size, costs, model, normalised, code_cost):
    filename = ""
    title = ""
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

    if normalised == False:
        title = "Number of occurrences of randomly simulated code costs per bracket for models generated using the " + model + " method\n Sample size = " + str(sample_size)
        filename = "_samples_" + model
    if normalised == True:
        title = "Number of occurrences of randomly simulated code costs per bracket for NORMALISED models generated using the " + model + " method\n Sample size = " + str(sample_size)
        filename = "_samples_" + model + "_NORM"    
 
    # create histogram
    line_label = model + " standard code cost"
    plt.figure(figsize=(15, 10))
    plt.hist(costs, bins=num_bins, edgecolor='white', linewidth=0.3)
    plt.axvline(x = code_cost, color = 'red', linestyle = '--', label = line_label)
    plt.text(code_cost + (0.01*(maxi-mini)), 2, rotation='vertical', s=line_label)
    plt.xlabel("Code cost")
    plt.ylabel("Number of occurrences")
    plt.title(title)
    plt.savefig("plots/hist" + f'{filename}.png')

    return 0

#make csv files of all matrices
def store_cost_matrices():
    filename = "matrices/SeqPredNN_cost_matrix.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in SeqPredNN_codon_matrix:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "matrices/Koonin_cost_matrix.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in Koonin_codon_matrix:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "matrices/SeqPredNN_cost_matrix_NORM.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in SeqPredNN_codon_matrix_NORM:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "matrices/Koonin_cost_matrix_NORM.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in Koonin_codon_matrix_NORM:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "matrices/Neutral_subst_cost_matrix.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in neutral_substitution_matrix:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1
    
    filename = "matrices/Higgs_cost_matrix.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in Higgs_codon_matrix:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

    filename = "matrices/Higgs_cost_matrix_NORM.csv"
    codon_string = " ," + ",".join(codons_excl_stop) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in Higgs_codon_matrix_NORM:
            file.write(codons_excl_stop[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

#start cost calculations
def calculate_SeqPredNN_Koonin_and_Higgs_matrices():
    for row in range(61):
        for cell in range(61):
            SeqPredNN_codon_matrix[row][cell] = get_cost(row, cell, "SeqPredNN", codons_per_aa, codons_excl_stop, False)
            Koonin_codon_matrix[row][cell] = get_cost(row, cell, "Koonin", codons_per_aa, codons_excl_stop, False)
            neutral_substitution_matrix[row][cell] = get_cost(row, cell, "neutral", codons_per_aa, codons_excl_stop, False)
            Higgs_codon_matrix[row][cell] = get_cost(row, cell, "Higgs", codons_per_aa, codons_excl_stop, False)

            if not isinstance(SeqPredNN_codon_matrix[row][cell], str):
                SeqPredNN_check.append(SeqPredNN_codon_matrix[row][cell])

            if not isinstance(Koonin_codon_matrix[row][cell], str):
                Koonin_check.append(Koonin_codon_matrix[row][cell])

            if not isinstance(neutral_substitution_matrix[row][cell], str):
                neutral_check.append(neutral_substitution_matrix[row][cell])

            if not isinstance(Higgs_codon_matrix[row][cell], str):
                Higgs_check.append(Higgs_codon_matrix[row][cell])

            SeqPredNN_codon_matrix_PLOTS[row][cell] = get_cost(row, cell, "SeqPredNN", codons_per_aa, codons_excl_stop, True)
            Koonin_codon_matrix_PLOTS[row][cell] = get_cost(row, cell, "Koonin", codons_per_aa, codons_excl_stop, True)
            neutral_substitution_matrix_PLOTS[row][cell] = get_cost(row, cell, "neutral", codons_per_aa, codons_excl_stop, True)
            Higgs_codon_matrix_PLOTS[row][cell] = get_cost(row, cell, "Higgs", codons_per_aa, codons_excl_stop, True)

def normalise_matrix(minimum, maximum, original_matrix, plot_matrix, norm_check):
    normalised = [[0 for x in range(61)] for y in range(61)]

    for row in range(len(original_matrix)):
        for cell in range(len(original_matrix)):
            if not isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
                plot_matrix[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
                norm_check.append(normalised[row][cell])
            if isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = original_matrix[row][cell]
                plot_matrix[row][cell] = -10

    return(normalised)

def normalise_matrix_SAMPLES(minimum, maximum, original_matrix):
    normalised = [[0 for x in range(61)] for y in range(61)]

    for row in range(len(original_matrix)):
        for cell in range(len(original_matrix)):
            if not isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
            if isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = original_matrix[row][cell]

    return(normalised)

def confusion_matrix(plot_matrix, original_matrix, minimum, maximum, filename, title):

    fig, ax = plt.subplots(figsize=(45, 38))
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.titlesize'] = 60

    cmap = matplotlib.colormaps.get_cmap('viridis')
    cmap.set_under(color='grey')
    img = ax.imshow(plot_matrix, cmap=cmap, vmin=minimum, vmax=maximum)

    plt.title(title, pad=30)
    ax.set_xticks(np.arange(61))
    ax.set_yticks(np.arange(61))
    ax.set_xticklabels(codons_excl_stop, rotation = 90)
    ax.set_yticklabels(codons_excl_stop)
    ax.tick_params(labelsize=30, pad=14, length=14, width=3)
    ax.tick_params(bottom=False, top=True, left=True, right=False)
    ax.tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=False)
    ax.set_xlabel('Mutant codon', fontsize=70, labelpad=24)
    ax.set_ylabel('Wild type codon', fontsize=60, labelpad=24)

    for i in range(61):
        for j in range(61):
            if not isinstance(original_matrix[i][j], str):
                ax.text(j, i, round(original_matrix[i][j], 3), ha="center", va="center", color="black")
            else:
                ax.text(j, i, original_matrix[i][j], ha="center", va="center", color="black")

    cbar = fig.colorbar(img, ax=ax, pad=0.01)
    cbar.ax.tick_params(labelsize=63, pad=14, length=14, width=3)
    fig.tight_layout()
    plt.savefig("plots/" + f'{filename}.png')
    plt.close()

def spearmans_rank_correlation_tests():

    filename = "stats/Spearmans_rank_correlation_tests.csv"
    with open(filename, mode="w") as file:
        file.write("code cost test, Null hypothesis, Alternative hypothesis, correlation, p-value\n")

        file.write("SeqPredNN vs Koonin, corr == 0, corr != 0," + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='two-sided').correlation) + ","  + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='two-sided').pvalue)+ "\n")
        file.write("SeqPredNN vs Koonin, corr >= 0, corr < 0," + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='less').correlation) + ","  + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='less').pvalue)+ "\n")
        file.write("SeqPredNN vs Koonin, corr <= 0, corr > 0," + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='greater').correlation) + ","  + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='greater').pvalue)+ "\n")

        file.write("SeqPredNN vs Higgs, corr == 0, corr != 0," + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='two-sided').correlation) + ","  + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='two-sided').pvalue)+ "\n")
        file.write("SeqPredNN vs Higgs, corr >= 0, corr < 0," + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='less').correlation) + ","  + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='less').pvalue)+ "\n")
        file.write("SeqPredNN vs Higgs, corr <= 0, corr > 0," + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='greater').correlation) + ","  + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='greater').pvalue)+ "\n")

        file.write("Higgs vs Koonin, corr == 0, corr != 0," + str(spearmanr(Koonin_check, Higgs_check, alternative='two-sided').correlation) + ","  + str(spearmanr(Koonin_check, Higgs_check, alternative='two-sided').pvalue)+ "\n")
        file.write("Higgs vs Koonin, corr >= 0, corr < 0," + str(spearmanr(Koonin_check, Higgs_check, alternative='less').correlation) + ","  + str(spearmanr(Koonin_check, Higgs_check, alternative='less').pvalue)+ "\n")
        file.write("Higgs vs Koonin, corr <= 0, corr > 0," + str(spearmanr(Koonin_check, Higgs_check, alternative='greater').correlation) + ","  + str(spearmanr(Koonin_check, Higgs_check, alternative='greater').pvalue)+ "\n")

#generate stats and store in csvs
def stats():
    filename = "stats/code_cost_stats.csv"
    with open(filename, mode="w") as file:
        file.write(" , raw code cost, min, max, NORM code cost, min NORM, max NORM\n")

        file.write("SeqPredNN," + str(SeqPredNN_code_cost) + "," + str(min(SeqPredNN_check)) + "," + str(max(SeqPredNN_check)) + "," + str(SeqPredNN_code_cost_NORM) + "," + str(min(SeqPredNN_check_NORM)) + "," + str(max(SeqPredNN_check_NORM)) + "\n")

        file.write("Koonin," + str(Koonin_code_cost) + "," + str(min(Koonin_check)) + "," + str(max(Koonin_check)) + "," + str(Koonin_code_cost_NORM) + "," + str(min(Koonin_check_NORM)) + "," + str(max(Koonin_check_NORM)) + "\n")

        file.write("Higgs," + str(Higgs_code_cost) + "," + str(min(Higgs_check)) + "," + str(max(Higgs_check)) + "," + str(Higgs_code_cost_NORM) + "," + str(min(Higgs_check_NORM)) + "," + str(max(Higgs_check_NORM)) + "\n")
        
    filename = "stats/Kolmogorov_Smirnov_tests.csv"
    with open(filename, mode="w") as file:
        file.write("sample code costs test, statistic, p-value\n")

        file.write("SeqPredNN vs Koonin," + str(ks_2samp(SeqPredNN_sample_code_costs, Koonin_sample_code_costs).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs, Koonin_sample_code_costs).pvalue)+ "\n")

        file.write("SeqPredNN NORM vs Koonin NORM," + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).pvalue)+ "\n")

        file.write("SeqPredNN vs Higgs," + str(ks_2samp(SeqPredNN_sample_code_costs, Higgs_sample_code_costs).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs, Higgs_sample_code_costs).pvalue)+ "\n")

        file.write("SeqPredNN NORM vs Higgs NORM," + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Higgs_sample_code_costs_NORM).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Higgs_sample_code_costs_NORM).pvalue)+ "\n")

        file.write("Higgs vs Koonin," + str(ks_2samp(Higgs_sample_code_costs, Koonin_sample_code_costs).statistic) + ","  + str(ks_2samp(Higgs_sample_code_costs, Koonin_sample_code_costs).pvalue)+ "\n")

        file.write("Higgs NORM vs Koonin NORM," + str(ks_2samp(Higgs_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).statistic) + ","  + str(ks_2samp(Higgs_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).pvalue)+ "\n")

    SeqPredNN_check.sort()
    Koonin_check.sort()
    Higgs_check.sort()

    #Spearman's rank correlation coefficient:
    spearmans_rank_correlation_tests()

#generate cost of primordial code as presented by Koonin and Novozhilov Fig 12a
def primordial_code_simulation(mode, output_matrix, matrix_check, output_matrix_NORM):
    #_______________________________________
    #       |   U   |   C   |   A   |   G   |
    #-------|-------|-------|-------|-------|
    #   U   |  Leu  |  Ser  |   *   |   *   |
    #-------|-------|-------|-------|-------|
    #   C   |  Leu  |  Pro  |  His  |  Arg  |
    #-------|-------|-------|-------|-------|
    #   A   |  Ile  |  Thr  |  Asn  |  Ser  |
    #-------|-------|-------|-------|-------|
    #   G   |  Val  |  Ala  |  Asp  |  Gly  |
    #_______|_______|_______|_______|_______|

    primordial_codons_per_aa = {"HIS": ['CA'], "ARG": ['CG'], "ASP": ['GA'], "ASN": ['AA'], "GLY": ['GG'], "ALA": ['GC'], "SER": ['UC', 'AG'], "THR": ['AC'], "PRO": ['CC'], "VAL": ['GU'], "ILE": ['AU'], "LEU": ['CU', 'UU']}
    primordial_codons = ['UU', 'UC', 'UA', 'UG', 'CU', 'CC', 'CA', 'CG', 'AU', 'AC', 'AA', 'AG', 'GU', 'GC', 'GA', 'GG']

    for row in range(16):
        for cell in range(16):
            true_codon = primordial_codons[row]
            mutant_codon = primordial_codons[cell]

            if true_codon != mutant_codon and not (true_codon[0:1] != mutant_codon[0:1] and true_codon[1:2] != mutant_codon[1:2]): #not the same codon and only differ by 1 position
                weight = 0
                aa_diff = 0
                # true_aa = ""
                # mutant_aa = ""
                if true_codon[0:1] != mutant_codon[0:1]: #differ at 1st position
                    if tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "ts":
                        weight = FreelandHurst_mutation_weights.get("1ts")
                    elif tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "tv":
                        weight = FreelandHurst_mutation_weights.get("1tv")
                if true_codon[1:2] != mutant_codon[1:2]: #differ at 2nd position
                    if tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "ts":
                        weight = FreelandHurst_mutation_weights.get("2ts")
                    elif tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "tv":
                        weight = FreelandHurst_mutation_weights.get("2tv")
                print("weight: " + str(weight))
                for key, value in primordial_codons_per_aa.items():
                    if true_codon in value:
                        true_aa = key
                    if mutant_codon in value:
                        mutant_aa = key

                
                # print("true: " + true_codon + ", " + true_aa)
                # print("mutant: " + mutant_codon + ", " + mutant_aa)
                
                if mode != "Higgs":
                    true_aa_index = amino_acids.index(true_aa)
                    mutant_aa_index = amino_acids.index(mutant_aa)
                elif mode == "Higgs":
                    true_aa_index = Higgs_amino_acid_order.index(true_aa)
                    mutant_aa_index = Higgs_amino_acid_order.index(mutant_aa)

                if mode == "SeqPredNN":
                    aa_diff = 1 - float(normalised_substitution_matrix[true_aa_index][mutant_aa_index])
                elif mode == "Koonin":
                    aa_diff = pow((polar_requirement_scale.get(true_aa) - polar_requirement_scale.get(mutant_aa)), 2)
                elif mode == "Higgs":
                    aa_diff = Higgs_distance_matrix[true_aa_index][mutant_aa_index]
                elif mode == "neutral":
                    aa_diff = 1

                print("aa diff: " + str(aa_diff))
                output_matrix[row][cell] = float(weight) * float(aa_diff)
                matrix_check.append(output_matrix[row][cell])
                print("cost: " + str(output_matrix[row][cell]))
                
            else:
                output_matrix[row][cell] = '-'
        
        print(output_matrix[row])
    
    print("NORM:")
    for row in range(16):
        for cell in range(16):
            if not isinstance(output_matrix[row][cell], str):
                output_matrix_NORM[row][cell] = normalise_cost(min(matrix_check), max(matrix_check), output_matrix[row][cell])
            else: 
                output_matrix_NORM[row][cell] = output_matrix[row][cell]
        print(output_matrix_NORM[row])
    
    filename = "matrices/" + mode + "_primordial_cost_matrix.csv"
    codon_string = " ," + ",".join(primordial_codons) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in output_matrix:
            file.write(primordial_codons[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1
    
    filename = "matrices/" + mode + "_primordial_cost_matrix_NORM.csv"
    codon_string = " ," + ",".join(primordial_codons) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in output_matrix_NORM:
            file.write(primordial_codons[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

#run calculations for SeqPredNN and Koonin matrices
calculate_SeqPredNN_Koonin_and_Higgs_matrices()

#create normalised codon mutation cost matrices
SeqPredNN_codon_matrix_NORM = normalise_matrix(min(SeqPredNN_check), max(SeqPredNN_check), SeqPredNN_codon_matrix, SeqPredNN_codon_matrix_NORM_PLOTS, SeqPredNN_check_NORM)
Koonin_codon_matrix_NORM = normalise_matrix(min(Koonin_check), max(Koonin_check), Koonin_codon_matrix, Koonin_codon_matrix_NORM_PLOTS, Koonin_check_NORM)
Higgs_codon_matrix_NORM = normalise_matrix(min(Higgs_check), max(Higgs_check), Higgs_codon_matrix, Higgs_codon_matrix_NORM_PLOTS, Higgs_check_NORM)

#get overall code cost for SeqPredNN and Koonin matrices
SeqPredNN_code_cost = get_code_cost(SeqPredNN_codon_matrix)
Koonin_code_cost = get_code_cost(Koonin_codon_matrix)
Higgs_code_cost = get_code_cost(Higgs_codon_matrix)

SeqPredNN_code_cost_NORM = get_code_cost(SeqPredNN_codon_matrix_NORM)
Koonin_code_cost_NORM = get_code_cost(Koonin_codon_matrix_NORM)
Higgs_code_cost_NORM = get_code_cost(Higgs_codon_matrix_NORM)

#generate random codon assignments for SeqPredNN model and calculate code costs for each assignment
generate_sample_set(10000, SeqPredNN_sample_code_costs, SeqPredNN_sample_code_costs_NORM, SeqPredNN_code_cost, SeqPredNN_code_cost_NORM, "SeqPredNN")
generate_sample_set(10000, Koonin_sample_code_costs, Koonin_sample_code_costs_NORM, Koonin_code_cost, Koonin_code_cost_NORM, "Koonin")
generate_sample_set(10000, Higgs_sample_code_costs, Higgs_sample_code_costs_NORM, Higgs_code_cost, Higgs_code_cost_NORM, "Higgs")

#generate .csv files
store_cost_matrices()

#generate stats and store in csvs
stats()


#generate plots
# confusion_matrix(SeqPredNN_codon_matrix_PLOTS, SeqPredNN_codon_matrix, min(SeqPredNN_check), max(SeqPredNN_check), "SeqPredNN_codon_matrix", "Confusion matrix of codon mutation costs calculated\nusing SeqPredNN amino acid frequencies")
# confusion_matrix(SeqPredNN_codon_matrix_NORM_PLOTS, SeqPredNN_codon_matrix_NORM, min(SeqPredNN_check_NORM), max(SeqPredNN_check_NORM), "SeqPredNN_codon_matrix_NORM", "Confusion matrix of normalised codon mutation costs calculated\nusing SeqPredNN amino acid frequencies")
# confusion_matrix(Koonin_codon_matrix_PLOTS, Koonin_codon_matrix, min(Koonin_check), max(Koonin_check), "Koonin_codon_matrix", "Confusion matrix of codon mutation costs calculated\nusing the Polarity Requirement Index")
# confusion_matrix(Koonin_codon_matrix_NORM_PLOTS, Koonin_codon_matrix_NORM, min(Koonin_check_NORM), max(Koonin_check_NORM), "Koonin_codon_matrix_NORM", "Confusion matrix of normalised codon mutation costs calculated\nusing the Polarity Requirement Index")
# confusion_matrix(neutral_substitution_matrix_PLOTS, neutral_substitution_matrix, min(neutral_check), max(neutral_check), "Neutral_subst_codon_matrix", "Confusion matrix of codon mutation costs calculated\nusing amino acid differences of 1 for all mutations")

#primordial simulations
# primordial_code_simulation("SeqPredNN", SeqPredNN_primordial_codon_matrix, SeqPredNN_primordial_check, SeqPredNN_primordial_codon_matrix_NORM)
# primordial_code_simulation("Koonin", Koonin_primordial_codon_matrix, Koonin_primordial_check, Koonin_primordial_codon_matrix_NORM)
# primordial_code_simulation("Higgs", Higgs_primordial_codon_matrix, Higgs_primordial_check, Higgs_primordial_codon_matrix_NORM)