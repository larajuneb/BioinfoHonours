from cmath import sqrt
from operator import index
import random
import csv
from statistics import mean
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
SeqPredNN_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
SeqPredNN_check = []
SeqPredNN_check_NORM = []
SeqPredNN_amino_acid_check = []
SeqPredNN_amino_acid_check_NORM = []
SeqPredNN_sample_code_costs = []
SeqPredNN_sample_code_costs_NORM = []

#Koonin codon matrix
Koonin_codon_matrix = [[0 for x in range(61)] for y in range(61)]
Koonin_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
Koonin_check = []
Koonin_check_NORM = []
Koonin_sample_code_costs = []
Koonin_sample_code_costs_NORM = []

#Higgs codon matrix
Higgs_codon_matrix = [[0 for x in range(61)] for y in range(61)]
Higgs_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
Higgs_check = []
Higgs_check_NORM = []
Higgs_sample_code_costs = []
Higgs_sample_code_costs_NORM = []

#neutral
neutral_codon_matrix = [[0 for x in range(61)] for y in range(61)]
neutral_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]
neutral_check = []
neutral_check_NORM = []
neutral_sample_code_costs = []
neutral_sample_code_costs_NORM = []

#SeqPredNN primordial
SeqPredNN_primordial_codon_matrix = [[0 for x in range(14)] for y in range(14)]
SeqPredNN_primordial_codon_matrix_NORM = [[0 for x in range(14)] for y in range(14)]
SeqPredNN_primordial_check = []
SeqPredNN_primordial_check_NORM = []
SeqPredNN_primordial_sample_code_costs = []
SeqPredNN_primordial_sample_code_costs_NORM = []

#Koonin primordial
Koonin_primordial_codon_matrix = [[0 for x in range(14)] for y in range(14)]
Koonin_primordial_codon_matrix_NORM = [[0 for x in range(14)] for y in range(14)]
Koonin_primordial_check = []
Koonin_primordial_check_NORM = []
Koonin_primordial_sample_code_costs = []
Koonin_primordial_sample_code_costs_NORM = []

#Higgs primordial
Higgs_primordial_codon_matrix = [[0 for x in range(14)] for y in range(14)]
Higgs_primordial_codon_matrix_NORM = [[0 for x in range(14)] for y in range(14)]
Higgs_primordial_check = []
Higgs_primordial_check_NORM = []
Higgs_primordial_sample_code_costs = []
Higgs_primordial_sample_code_costs_NORM = []

#neutral primordial
neutral_primordial_codon_matrix = [[0 for x in range(14)] for y in range(14)]
neutral_primordial_codon_matrix_NORM = [[0 for x in range(14)] for y in range(14)]
neutral_primordial_check = []
neutral_primordial_check_NORM = []
neutral_primordial_sample_code_costs = []
neutral_primordial_sample_code_costs_NORM = []

#amino acid primordial
amino_acid_primordial_codon_matrix = [[0 for x in range(14)] for y in range(14)]
amino_acid_primordial_codon_matrix_NORM = [[0 for x in range(14)] for y in range(14)]
amino_acid_primordial_check = []
amino_acid_primordial_check_NORM = []
amino_acid_primordial_sample_code_costs = []
amino_acid_primordial_sample_code_costs_NORM = []

#SeqPredNN amino acid substitutions only
amino_acid_codon_matrix = [[0 for x in range(61)] for y in range(61)]
amino_acid_codon_matrix_NORM = [[0 for x in range(61)] for y in range(61)]

nucleotides = ['U', 'C', 'A', 'G']
codons = []
codons_excl_stop = []
codons_leftover = []
codon = ""
stop = ['UAA', 'UAG', 'UGA']
disregard = ['UA', 'UG']

#generate all codons
for first in nucleotides:
    for second in nucleotides:
        for third in nucleotides:
            codon = first + second + third
            codons.append(codon)
            codons_leftover.append(codon)
            if codon not in stop:
                codons_excl_stop.append(codon)

number_of_codons_per_aa = {"HIS": 2, "ARG": 6, "LYS": 2, "GLN": 2, "GLU": 2, "ASP": 2, "ASN": 2, "GLY": 4, "ALA": 4, "SER": 6, "THR": 4, "PRO": 4, "CYS": 2, "VAL": 4, "ILE": 3, "MET": 1, "LEU": 6, "PHE": 2, "TYR": 2, "TRP": 1, "stop": 3}
amino_acids = ["HIS", "ARG", "LYS", "GLN", "GLU", "ASP", "ASN", "GLY", "ALA", "SER", "THR", "PRO", "CYS", "VAL", "ILE", "MET", "LEU", "PHE", "TYR", "TRP"]
codons_per_aa = {"HIS": ['CAU', 'CAC'], "ARG": ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], "LYS": ['AAA', 'AAG'], "GLN": ['CAA', 'CAG'], "GLU": ['GAA', 'GAG'], "ASP": ['GAU', 'GAC'], "ASN": ['AAU', 'AAC'], "GLY": ['GGU', 'GGC', 'GGA', 'GGG'], "ALA": ['GCU', 'GCC', 'GCA', 'GCG'], "SER": ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], "THR": ['ACU', 'ACC', 'ACA', 'ACG'], "PRO": ['CCU', 'CCC', 'CCA', 'CCG'], "CYS": ['UGU', 'UGC'], "VAL": ['GUU', 'GUC', 'GUA', 'GUG'], "ILE": ['AUU', 'AUC', 'AUA'], "MET": ['AUG'], "LEU": ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'], "PHE": ['UUU', 'UUC'], "TYR": ['UAU', 'UAC'], "TRP": ['UGG'], "stop": ['UAA', 'UAG', 'UGA']}

Higgs_amino_acid_order = ['PHE','LEU','ILE','MET','VAL','SER','PRO','THR','ALA','TYR','HIS','GLN','ASN','LYS','ASP','GLU','CYS','TRP','ARG','GLY']


primordial_number_of_codons_per_aa = {"HIS": 1, "ARG": 1, "ASP": 1, "ASN": 1, "GLY": 1, "ALA": 1, "SER": 2, "THR": 1, "PRO": 1, "VAL": 1, "ILE": 1, "LEU": 2, "disregard": 2}
primordial_amino_acids = ["HIS", "ARG", "ASP", "ASN", "GLY", "ALA", "SER", "THR", "PRO", "VAL", "ILE", "LEU"]
primordial_codons_per_aa = {"HIS": ['CA'], "ARG": ['CG'], "ASP": ['GA'], "ASN": ['AA'], "GLY": ['GG'], "ALA": ['GC'], "SER": ['UC', 'AG'], "THR": ['AC'], "PRO": ['CC'], "VAL": ['GU'], "ILE": ['AU'], "LEU": ['CU', 'UU'], "disregard": ['UA', 'UG']}
primordial_codons = ['UU', 'UC', 'UA', 'UG', 'CU', 'CC', 'CA', 'CG', 'AU', 'AC', 'AA', 'AG', 'GU', 'GC', 'GA', 'GG'] #UA and UG should be excluded
primordial_codons_excl_stop = ['UU', 'UC', 'CU', 'CC', 'CA', 'CG', 'AU', 'AC', 'AA', 'AG', 'GU', 'GC', 'GA', 'GG']


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
def get_cost(true_codon_index, mutant_codon_index, model, codon_dict, codons_without_stop, plot):
    #get codon string name from it's index
    true_codon = codons_without_stop[true_codon_index]
    mutant_codon = codons_without_stop[mutant_codon_index]
    #get amino acid name from codon name
    true_aa = get_aa_for_codon(true_codon, codon_dict)
    mutant_aa = get_aa_for_codon(mutant_codon, codon_dict)
    #get amino acid name from amino acid index
    if model != "Higgs":
        true_aa_index = amino_acids.index(true_aa)
        mutant_aa_index = amino_acids.index(mutant_aa)
    elif model == "Higgs":
        true_aa_index = Higgs_amino_acid_order.index(true_aa)
        mutant_aa_index = Higgs_amino_acid_order.index(mutant_aa)

    aa_difference = 0
    if model == "SeqPredNN":
        aa_difference = 1 - float(normalised_substitution_matrix[true_aa_index][mutant_aa_index])
    elif model == "Koonin":
        aa_difference = pow((polar_requirement_scale.get(true_aa) - polar_requirement_scale.get(mutant_aa)), 2)
    elif model == "Higgs":
        aa_difference = Higgs_distance_matrix[true_aa_index][mutant_aa_index]
    elif model == "neutral":
        aa_difference = 1

    codon_mutation_prob = get_codon_mutation_prob(true_codon, mutant_codon)
    if codon_mutation_prob == '*':
        cost = '*'
    elif codon_mutation_prob == '-':
        cost = '-'
    elif model == "amino-acid":
        cost = 1 - float(normalised_substitution_matrix[true_aa_index][mutant_aa_index])
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
    for row in range(len(cost_array)):
        for cell in range(len(cost_array)):
            if not isinstance(cost_array[row][cell], str):
                code_cost += cost_array[row][cell]
    return code_cost

def regenerate_codon_matrix(array, mode):
    for first in nucleotides:
        for second in nucleotides:
            if mode == "standard":
                for third in nucleotides:
                    codon = first + second + third
                    array.append(codon)
            elif mode == "primordial":
                codon = first + second
                array.append(codon)
    return(array)

#check if a list of costs in normally distributed using the Kolmogorov-Smirnov Test
def is_it_gaussian(dataset):
    kstest(dataset, 'norm')
    return 0

#generate 10,000 random assignments of codons and calculate costs for each
def generate_sample_set(sample_size, sample_code_costs, sample_code_costs_NORM, mode_code_cost, mode_code_cost_NORM, model, neutral_cost, neutral_cost_NORM, mode, matrix_length):
    random_codon_assignments = {} #similar to codons_per_aa dict, but instead of true codon assignments, the codons are assigned randomly to amino acids
    leftover = []
    exclude = []
    code_costs = []
    norm_code_costs = []
    no_stop_codons = []
    if mode == "standard":
        number_of_codons_per_aa_SAMPLE = number_of_codons_per_aa
    elif mode == "primordial":
        number_of_codons_per_aa_SAMPLE = primordial_number_of_codons_per_aa

    for i in range(sample_size):
        temp_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        norm_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        matrix_min_max_check = []
        minimum = 0
        maximum = 0
        leftover = regenerate_codon_matrix(leftover, mode)
        no_stop_codons = regenerate_codon_matrix(no_stop_codons, mode)
        for key, value in number_of_codons_per_aa_SAMPLE.items():
            random_codon_assignments[key] = [] #a new random codon assignment for each sample, key = amino acid, value = array of codons
            for i in range(number_of_codons_per_aa_SAMPLE[key]): #loop through as many times as there are codons for that amino acid
                codon = random.choice(leftover) #choose a random codon for the amino acid 
                leftover.remove(codon) #remove that codon from the list of available codons
                random_codon_assignments[key].append(codon) #add the randomly chose codon to the list of codons for that amino acid
        if mode == "standard":
            exclude = random_codon_assignments.get("stop")
        elif mode == "primordial":
            exclude = random_codon_assignments.get("disregard")
        #random codon assignment have been made, now calculate cost matrix and overall cost
        for codon in exclude:
            no_stop_codons.remove(codon)
        
        for row in range(matrix_length):
            for cell in range(matrix_length):
                temp_cost_matrix[row][cell] = get_cost(row, cell, model, random_codon_assignments, no_stop_codons, False)

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
        temp_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        norm_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        
    sample_set_stats(code_costs, norm_code_costs, model)

    plot_samples(sample_size, code_costs, model, False, mode_code_cost, neutral_cost)
    plot_samples(sample_size, norm_code_costs, model, True, mode_code_cost_NORM, neutral_cost_NORM)

    code_costs.clear()
    norm_code_costs.clear()

def sample_set_stats(costs, norm_costs, mode):
    series = pd.Series(costs)
    series_norm = pd.Series(norm_costs)

    filename = "TRY/stats/" + mode + "_code_cost_samples_stats.csv"
    with open(filename, mode="w") as file:
        file.write("mode, count, mean, std, min, 25%, 50%, 75%, max\n")

        file.write(mode + "," + str(series.describe()[0]) + "," + str(series.describe()[1]) + "," + str(series.describe()[2]) + "," + str(series.describe()[3]) + "," + str(series.describe()[4]) + "," + str(series.describe()[5]) + "," + str(series.describe()[6]) + "," + str(series.describe()[7]) + "\n")

        file.write(mode + "NORM," + str(series_norm.describe()[0]) + "," + str(series_norm.describe()[1]) + "," + str(series_norm.describe()[2]) + "," + str(series_norm.describe()[3]) + "," + str(series_norm.describe()[4]) + "," + str(series_norm.describe()[5]) + "," + str(series_norm.describe()[6]) + "," + str(series_norm.describe()[7]) + "\n")

#plot bar graphs for samples produced
def plot_samples(sample_size, costs, model, normalised, code_cost, neutral_cost):
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
    plt.axvline(x = code_cost, color = 'red', linestyle = '--', label = "hi")
    plt.text(code_cost + (0.01*(maxi-mini)), 2, rotation='vertical', s=line_label)

    line_label = "Neutral substitutions standard code cost"
    plt.axvline(x = neutral_cost, color = 'green', linestyle = '--', label = line_label)
    plt.text(neutral_cost + (0.01*(maxi-mini)), 2, rotation='vertical', s=line_label)
    plt.xlabel("Code cost")
    plt.ylabel("Number of occurrences")
    plt.title(title)
    plt.savefig("TRY/plots/histograms/hist" + f'{filename}.png')

    return 0

#make csv files of all matrices
def store_cost_matrices(mode, matrix):
    filename = "TRY/matrices/" + mode + ".csv"
    if "primordial" in mode:
        codon_set = primordial_codons_excl_stop
    else:
        codon_set = codons_excl_stop
    codon_string = " ," + ",".join(codon_set) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in matrix:
            file.write(codon_set[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

#start cost calculations
def calculate_cost_matrix(codon_matrix, check, model, mode):
    for row in range(len(codon_matrix)):
        for cell in range(len(codon_matrix)):
            if mode == "standard":
                codon_matrix[row][cell] = get_cost(row, cell, model, codons_per_aa, codons_excl_stop, False)
            elif mode == "primordial":
                codon_matrix[row][cell] = get_cost(row, cell, model, primordial_codons_per_aa, primordial_codons_excl_stop, False)

            if not isinstance(codon_matrix[row][cell], str):
                check.append(codon_matrix[row][cell])

def normalise_matrix(minimum, maximum, original_matrix, norm_check):
    normalised = [[0 for x in range(len(original_matrix))] for y in range(len(original_matrix))]

    for row in range(len(original_matrix)):
        for cell in range(len(original_matrix)):
            if not isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
                # plot_matrix[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
                norm_check.append(normalised[row][cell])
            if isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = original_matrix[row][cell]
                # plot_matrix[row][cell] = -10

    return(normalised)

def normalise_matrix_SAMPLES(minimum, maximum, original_matrix):
    normalised = [[0 for x in range(len(original_matrix))] for y in range(len(original_matrix))]

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

    for i in range(len(original_matrix)):
        for j in range(len(original_matrix)):
            if not isinstance(original_matrix[i][j], str):
                ax.text(j, i, round(original_matrix[i][j], 3), ha="center", va="center", color="black")
            else:
                ax.text(j, i, original_matrix[i][j], ha="center", va="center", color="black")

    cbar = fig.colorbar(img, ax=ax, pad=0.01)
    cbar.ax.tick_params(labelsize=63, pad=14, length=14, width=3)
    fig.tight_layout()
    plt.savefig("TRY/plots/" + f'{filename}.png')
    plt.close()

def spearmans_rank_correlation_tests():

    filename = "TRY/stats/Spearmans_rank_correlation_tests.csv"
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
    filename = "TRY/stats/code_cost_stats.csv"
    with open(filename, mode="w") as file:
        file.write(" , raw code cost, mean, min, max, NORM code cost, mean NORM, min NORM, max NORM\n")

        file.write("SeqPredNN," + str(SeqPredNN_code_cost) + "," + str(mean(SeqPredNN_check)) + "," + str(min(SeqPredNN_check)) + "," + str(max(SeqPredNN_check)) + "," + str(SeqPredNN_code_cost_NORM) + "," + str(mean(SeqPredNN_check_NORM)) + "," + str(min(SeqPredNN_check_NORM)) + "," + str(max(SeqPredNN_check_NORM)) + "\n")

        file.write("Koonin," + str(Koonin_code_cost) + "," + str(mean(Koonin_check)) + "," + str(min(Koonin_check)) + "," + str(max(Koonin_check)) + "," + str(Koonin_code_cost_NORM) + "," + str(mean(Koonin_check_NORM)) + "," + str(min(Koonin_check_NORM)) + "," + str(max(Koonin_check_NORM)) + "\n")

        file.write("Higgs," + str(Higgs_code_cost) + "," + str(mean(Higgs_check)) + "," + str(min(Higgs_check)) + "," + str(max(Higgs_check)) + "," + str(Higgs_code_cost_NORM) + "," + str(mean(Higgs_check_NORM)) + "," + str(min(Higgs_check_NORM)) + "," + str(max(Higgs_check_NORM)) + "\n")
        
    filename = "TRY/stats/Kolmogorov_Smirnov_tests.csv"
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

#run calculations for SeqPredNN, Koonin, Higgs, neutral and amino-acid substitution matrices
calculate_cost_matrix(SeqPredNN_codon_matrix, SeqPredNN_check, "SeqPredNN", "standard")
calculate_cost_matrix(Koonin_codon_matrix, Koonin_check, "Koonin", "standard")
calculate_cost_matrix(Higgs_codon_matrix, Higgs_check, "Higgs", "standard")
calculate_cost_matrix(neutral_codon_matrix, neutral_check, "neutral", "standard")
calculate_cost_matrix(amino_acid_codon_matrix, SeqPredNN_amino_acid_check, "amino-acid", "standard")
print(SeqPredNN_amino_acid_check)

#create normalised codon mutation cost matrices
SeqPredNN_codon_matrix_NORM = normalise_matrix(min(SeqPredNN_check), max(SeqPredNN_check), SeqPredNN_codon_matrix, SeqPredNN_check_NORM)
Koonin_codon_matrix_NORM = normalise_matrix(min(Koonin_check), max(Koonin_check), Koonin_codon_matrix, Koonin_check_NORM)
Higgs_codon_matrix_NORM = normalise_matrix(min(Higgs_check), max(Higgs_check), Higgs_codon_matrix, Higgs_check_NORM)
neutral_codon_matrix_NORM = normalise_matrix(min(neutral_check), max(neutral_check), neutral_codon_matrix, neutral_check_NORM)
amino_acid_codon_matrix_NORM = normalise_matrix(min(SeqPredNN_amino_acid_check), max(SeqPredNN_amino_acid_check), amino_acid_codon_matrix, SeqPredNN_amino_acid_check_NORM)

#get overall code cost for all matrices
SeqPredNN_code_cost = get_code_cost(SeqPredNN_codon_matrix)
Koonin_code_cost = get_code_cost(Koonin_codon_matrix)
Higgs_code_cost = get_code_cost(Higgs_codon_matrix)
neutral_code_cost = get_code_cost(neutral_codon_matrix)
amino_acid_code_cost = get_code_cost(amino_acid_codon_matrix)

SeqPredNN_code_cost_NORM = get_code_cost(SeqPredNN_codon_matrix_NORM)
Koonin_code_cost_NORM = get_code_cost(Koonin_codon_matrix_NORM)
Higgs_code_cost_NORM = get_code_cost(Higgs_codon_matrix_NORM)
neutral_code_cost_NORM = get_code_cost(neutral_codon_matrix_NORM)
amino_acid_code_cost_NORM = get_code_cost(amino_acid_codon_matrix_NORM)

generate_sample_set(100, SeqPredNN_sample_code_costs, SeqPredNN_sample_code_costs_NORM, SeqPredNN_code_cost, SeqPredNN_code_cost_NORM, "SeqPredNN", neutral_code_cost, neutral_code_cost_NORM, "standard", 61)
generate_sample_set(100, Koonin_sample_code_costs, Koonin_sample_code_costs_NORM, Koonin_code_cost, Koonin_code_cost_NORM, "Koonin", neutral_code_cost, neutral_code_cost_NORM, "standard", 61)
generate_sample_set(100, Higgs_sample_code_costs, Higgs_sample_code_costs_NORM, Higgs_code_cost, Higgs_code_cost_NORM, "Higgs", neutral_code_cost, neutral_code_cost_NORM, "standard", 61)
generate_sample_set(100, neutral_sample_code_costs, neutral_sample_code_costs_NORM, neutral_code_cost, neutral_code_cost_NORM, "neutral", neutral_code_cost, neutral_code_cost_NORM, "standard", 61)

calculate_cost_matrix(SeqPredNN_primordial_codon_matrix, SeqPredNN_primordial_check, "SeqPredNN", "primordial")
calculate_cost_matrix(Koonin_primordial_codon_matrix, Koonin_primordial_check, "Koonin", "primordial")
calculate_cost_matrix(Higgs_primordial_codon_matrix, Higgs_primordial_check, "Higgs", "primordial")
calculate_cost_matrix(neutral_primordial_codon_matrix, neutral_primordial_check, "neutral", "primordial")
calculate_cost_matrix(amino_acid_primordial_codon_matrix, amino_acid_primordial_check, "amino-acid", "primordial")

#get overall code cost for all primordial matrices
SeqPredNN_primordial_code_cost = get_code_cost(SeqPredNN_primordial_codon_matrix)
Koonin_primordial_code_cost = get_code_cost(Koonin_primordial_codon_matrix)
Higgs_primordial_code_cost = get_code_cost(Higgs_primordial_codon_matrix)
neutral_primordial_code_cost = get_code_cost(neutral_primordial_codon_matrix)
amino_acid_primordial_code_cost = get_code_cost(amino_acid_codon_matrix)

#create normalised codon mutation cost matrices
SeqPredNN_primordial_codon_matrix_NORM = normalise_matrix(min(SeqPredNN_primordial_check), max(SeqPredNN_primordial_check), SeqPredNN_primordial_codon_matrix, SeqPredNN_primordial_check_NORM)
Koonin_primordial_codon_matrix_NORM = normalise_matrix(min(Koonin_primordial_check), max(Koonin_primordial_check), Koonin_primordial_codon_matrix, Koonin_primordial_check_NORM)
Higgs_primordial_codon_matrix_NORM = normalise_matrix(min(Higgs_primordial_check), max(Higgs_primordial_check), Higgs_primordial_codon_matrix, Higgs_primordial_check_NORM)
neutral_primordial_codon_matrix_NORM = normalise_matrix(min(neutral_primordial_check), max(neutral_primordial_check), neutral_primordial_codon_matrix, neutral_primordial_check_NORM)
amino_acid_primordial_codon_matrix_NORM = normalise_matrix(min(amino_acid_primordial_check), max(amino_acid_primordial_check), amino_acid_primordial_codon_matrix, amino_acid_primordial_check_NORM)

SeqPredNN_primordial_code_cost_NORM = get_code_cost(SeqPredNN_primordial_codon_matrix_NORM)
Koonin_primordial_code_cost_NORM = get_code_cost(Koonin_primordial_codon_matrix_NORM)
Higgs_primordial_code_cost_NORM = get_code_cost(Higgs_primordial_codon_matrix_NORM)
neutral_primordial_code_cost_NORM = get_code_cost(neutral_primordial_codon_matrix_NORM)
amino_acid_primordial_code_cost_NORM = get_code_cost(amino_acid_primordial_codon_matrix_NORM)

generate_sample_set(100, SeqPredNN_primordial_sample_code_costs, SeqPredNN_primordial_sample_code_costs_NORM, SeqPredNN_primordial_code_cost, SeqPredNN_primordial_code_cost_NORM, "SeqPredNN", neutral_primordial_code_cost, neutral_primordial_code_cost_NORM, "primordial", 14)
generate_sample_set(100, Koonin_primordial_sample_code_costs, Koonin_primordial_sample_code_costs_NORM, Koonin_primordial_code_cost, Koonin_primordial_code_cost_NORM, "Koonin", neutral_primordial_code_cost, neutral_primordial_code_cost_NORM, "primordial", 14)
generate_sample_set(100, Higgs_primordial_sample_code_costs, Higgs_primordial_sample_code_costs_NORM, Higgs_primordial_code_cost, Higgs_primordial_code_cost_NORM, "Higgs", neutral_primordial_code_cost, neutral_primordial_code_cost_NORM, "primordial", 14)
generate_sample_set(100, neutral_primordial_sample_code_costs, neutral_primordial_sample_code_costs_NORM, neutral_primordial_code_cost, neutral_primordial_code_cost_NORM, "neutral", neutral_primordial_code_cost, neutral_primordial_code_cost_NORM, "primordial", 14)
generate_sample_set(100, amino_acid_primordial_sample_code_costs, amino_acid_primordial_sample_code_costs_NORM, amino_acid_primordial_code_cost, amino_acid_primordial_code_cost_NORM, "amino-acid", neutral_primordial_code_cost, neutral_primordial_code_cost_NORM, "primordial", 14)

#generate .csv files
store_cost_matrices("SeqPredNN_cost_matrix", SeqPredNN_codon_matrix)
store_cost_matrices("SeqPredNN_cost_matrix_NORM", SeqPredNN_codon_matrix_NORM)
store_cost_matrices("Koonin_cost_matrix", Koonin_codon_matrix)
store_cost_matrices("Koonin_cost_matrix_NORM", Koonin_codon_matrix_NORM)
store_cost_matrices("Higgs_cost_matrix", Higgs_codon_matrix)
store_cost_matrices("Higgs_cost_matrix_NORM", Higgs_codon_matrix_NORM)
store_cost_matrices("Neutral_subst_cost_matrix", neutral_codon_matrix)
store_cost_matrices("Neutral_subst_cost_matrix_NORM", neutral_codon_matrix_NORM)
store_cost_matrices("Amino_acid_codon_matrix", amino_acid_codon_matrix)
store_cost_matrices("Amino_acid_codon_matrix_NORM", amino_acid_codon_matrix_NORM)


store_cost_matrices("primordial/SeqPredNN_primordial_cost_matrix", SeqPredNN_primordial_codon_matrix)
store_cost_matrices("primordial/SeqPredNN_primordial_cost_matrix_NORM", SeqPredNN_primordial_codon_matrix_NORM)
store_cost_matrices("primordial/Koonin_primordial_cost_matrix", Koonin_primordial_codon_matrix)
store_cost_matrices("primordial/Koonin_primordial_cost_matrix_NORM", Koonin_primordial_codon_matrix_NORM)
store_cost_matrices("primordial/Higgs_primordial_cost_matrix", Higgs_primordial_codon_matrix)
store_cost_matrices("primordial/Higgs_primordial_cost_matrix_NORM", Higgs_primordial_codon_matrix_NORM)
store_cost_matrices("primordial/Neutral_primordial_subst_cost_matrix", neutral_primordial_codon_matrix)
store_cost_matrices("primordial/Neutral_primordial_subst_cost_matrix_NORM", neutral_primordial_codon_matrix_NORM)
store_cost_matrices("primordial/Amino_primordial_acid_codon_matrix", amino_acid_primordial_codon_matrix)
store_cost_matrices("primordial/Amino_primordial_acid_codon_matrix_NORM", amino_acid_primordial_codon_matrix_NORM)

#generate stats and store in csvs
# stats()