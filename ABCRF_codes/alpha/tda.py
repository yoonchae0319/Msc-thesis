import numpy as np
import sys
import csv
from ripser import Rips
import subprocess

import scipy
import sgtl
import sgtl.graph
import sgtl.spectrum
import sgtl.random
import sgtl.clustering
import scipy.sparse as sp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sns
import matplotlib.patches as patches
from pydtmc import MarkovChain


#transition_time = float(sys.argv[1])
alpha = float(sys.argv[1])           # evolutionary parameter from SIMgen.py
popsize = float(sys.argv[2])
sim_number = int(sys.argv[3])    # simulation index from loop

rips = Rips(maxdim=2)

# Load all distance matrices for this simulation
distance_matrices = []
j = 1
while True:
    try:
        filename = f"distanceMatrix_{sim_number}_{j}.txt"
        dist_matrix = np.loadtxt(filename)
        distance_matrices.append(dist_matrix)
        j += 1
    except OSError:
        break  # no more files for this simulation


##obtaining psi, the mean length of the zero homology group barcodes
def Psi(window, distanceMatrices):
    diagrams = rips.fit_transform(distanceMatrices[window], distance_matrix=True)
    h0 = diagrams[0]
    lengths = h0[:, 1] - h0[:, 0]  # death - birth
    finite_lengths = lengths[:-1]  # Exclude infinite bar
    if len(finite_lengths) == 0:
        return 0.0
    else:
        return np.mean(finite_lengths)

def Psi2(window, distanceMatrices):
    diagrams = rips.fit_transform(distanceMatrices[window], distance_matrix=True)
    h1 = diagrams[1]
    lengths = h1[:, 1] - h1[:, 0]  # death - birth
    if len(lengths) == 0:
        return 0.0
    else:
        return np.mean(lengths)

def Psi3(window, distanceMatrices):
    diagrams = rips.fit_transform(distanceMatrices[window], distance_matrix=True)
    h2 = diagrams[2]
    lengths = h2[:, 1] - h2[:, 0]  # death - birth
    if len(lengths) == 0:
        return 0.0
    else:
        return np.mean(lengths)

def b1(window, distanceMatrices):
    diagrams = rips.fit_transform(distanceMatrices[window], distance_matrix=True)
    h1 = diagrams[1]
    return (len(h1))

def b2(window, distanceMatrices):
    diagrams = rips.fit_transform(distanceMatrices[window], distance_matrix=True)
    h2 = diagrams[2]
    return (len(h2))

def tpi(window, distanceMatrices, n):
    return (np.sum(distanceMatrices[window]) / (n * (n - 1)))

def tajimas_D(pi, S, n):
    if S == 0:
        return 0.0  # No segregating sites, D is 0 or undefined
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i ** 2) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n ** 2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / a1 ** 2)
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    var = e1 * S + e2 * S * (S - 1)
    if var == 0:
        return 0.0  # Avoid divide by zero
    theta_w = S / a1
    D = (pi - theta_w) / np.sqrt(var)
    return D

def cluster(window, distance_matrices, num_clusters):
    cluster_ = []
    D_matrix = distance_matrices[window]
    sigma = np.mean(D_matrix)
    if sigma == 0 or np.isnan(sigma):
        print(f"Warning: sigma=0 or NaN at window {window}")
        return 0  # Or a safe fallback like 'num_cluster'
    A = np.exp(-D_matrix**2 / (2 * sigma**2))
    np.fill_diagonal(A, 0)  # no self-loops

    #Convert your dense NumPy array to a sparse CSR matrix
    D_sparse = sp.csr_matrix(A)

    # Now create the graph
    graph = sgtl.graph.Graph(D_sparse)
    clusters = sgtl.clustering.spectral_clustering(graph, num_clusters=num_clusters)
    for i in range(0,len(clusters)):
        cluster_.append(len(clusters[i]))
    return np.max(cluster_)

def clusterv(window, distance_matrices, num_clusters):
    clusterv_ = []
    D_matrix = distance_matrices[window]
    sigma = np.mean(D_matrix)
    if sigma == 0 or np.isnan(sigma):
        print(f"Warning: sigma=0 or NaN at window {window}")
        return 0  # Or a safe fallback like 'num_cluster'
    A = np.exp(-D_matrix**2 / (2 * sigma**2))
    np.fill_diagonal(A, 0)  # no self-loops

    #Convert your dense NumPy array to a sparse CSR matrix
    D_sparse = sp.csr_matrix(A)

    # Now create the graph
    graph = sgtl.graph.Graph(D_sparse)
    clusters = sgtl.clustering.spectral_clustering(graph, num_clusters=num_clusters)
    for i in range(0,len(clusters)):
        clusterv_.append(len(clusters[i]))
    return np.var(clusterv_)


def compute_transition_difference(data_array, log_edges):
    max_val = np.max(data_array)
    if max_val == 0:
        return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

    log_edges_norm = (log_edges - log_edges.min()) / (log_edges.max() - log_edges.min())
    scaled_edges = log_edges_norm * max_val

    num_states = len(scaled_edges) - 1
    transition_matrix = np.zeros((num_states, num_states))

    if len(data_array.shape) == 1:
        series = data_array
        states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
        states = pd.Series(states).fillna(0).astype(int)
        for i in range(len(states) - 1):
            current = states[i]
            next_ = states[i + 1]
            transition_matrix[current, next_] += 1
    else:
        for col in range(data_array.shape[1]):
            series = data_array[:, col]
            states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
            states = pd.Series(states).fillna(0).astype(int)
            for i in range(len(states) - 1):
                current = states[i]
                next_ = states[i + 1]
                transition_matrix[current, next_] += 1

    row_sums = transition_matrix.sum(axis=1, keepdims=True)
    transition_probs = np.nan_to_num(transition_matrix / row_sums)
    mat = np.array(transition_probs)

    a = num_states/5
    mat_sum = []
    for i in range(0,5):
        for ii in range(0,5):
            mat_sum.append(mat[int(i*a):int((i+1)*a), int(ii*a):int((ii+1)*a)].sum())

    var_lastcolumn = mat[:, :-1].var()
    var_lastrow = mat[:-1, :].var()


    eigenvalues, eigenvectors = np.linalg.eig(mat)
    eigen_var = eigenvalues.var()
    zero_count = np.count_nonzero(mat == 0)
    return mat_sum + [var_lastcolumn, var_lastrow, eigen_var, zero_count]

def compute_transition_difference2(data_array, log_edges): #for deta values and relatedness measures
    max_val = np.max(data_array)
    if max_val == 0:
        return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

    log_edges_norm = (log_edges.max()-log_edges[::-1]) / (log_edges.max() - log_edges.min())
    scaled_edges = log_edges_norm * max_val

    num_states = len(scaled_edges) - 1
    transition_matrix = np.zeros((num_states, num_states))

    if len(data_array.shape) == 1:
        series = data_array
        states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
        states = pd.Series(states).fillna(0).astype(int)
        for i in range(len(states) - 1):
            current = states[i]
            next_ = states[i + 1]
            transition_matrix[current, next_] += 1
    else:
        for col in range(data_array.shape[1]):
            series = data_array[:, col]
            states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
            states = pd.Series(states).fillna(0).astype(int)
            for i in range(len(states) - 1):
                current = states[i]
                next_ = states[i + 1]
                transition_matrix[current, next_] += 1

    row_sums = transition_matrix.sum(axis=1, keepdims=True)
    transition_probs = np.nan_to_num(transition_matrix / row_sums)
    mat = np.array(transition_probs)

    a = num_states/5
    mat_sum = []
    for i in range(0,5):
        for ii in range(0,5):
            mat_sum.append(mat[int(i*a):int((i+1)*a), int(ii*a):int((ii+1)*a)].sum())

    var_lastcolumn = mat[:, :-1].var()
    var_lastrow = mat[:-1, :].var()


    eigenvalues, eigenvectors = np.linalg.eig(mat)
    eigen_var = eigenvalues.var()
    zero_count = np.count_nonzero(mat == 0)
    return mat_sum + [var_lastcolumn, var_lastrow, eigen_var, zero_count]


def compute_ht(data_array, log_edges):
    max_val = np.max(data_array)
    if max_val == 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    log_edges_norm = (log_edges - log_edges.min()) / (log_edges.max() - log_edges.min())
    scaled_edges = log_edges_norm * max_val

    if len(np.unique(scaled_edges)) != len(scaled_edges):
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    num_states = len(scaled_edges) - 1
    transition_matrix = np.zeros((num_states, num_states))

    if len(data_array.shape) == 1:
        series = data_array
        states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
        states = pd.Series(states).fillna(0).astype(int)
        for i in range(len(states) - 1):
            current = states[i]
            next_ = states[i + 1]
            transition_matrix[current, next_] += 1
    else:
        for col in range(data_array.shape[1]):
            series = data_array[:, col]
            states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
            states = pd.Series(states).fillna(0).astype(int)
            for i in range(len(states) - 1):
                current = states[i]
                next_ = states[i + 1]
                transition_matrix[current, next_] += 1

    row_sums = transition_matrix.sum(axis=1, keepdims=True)
    transition_probs = np.nan_to_num(transition_matrix / row_sums)
    mat = np.array(transition_probs)
    for i in range(mat.shape[0]):
        if mat[i].sum() == 0:
            mat[i] = np.full(mat.shape[1], 1.0 / mat.shape[1])

    mc = MarkovChain(mat, ['1', '2', '3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20',
                           '21', '22', '23', '24', '25', '26', '27', '28', '29','30', '31', '32', '33', '34', '35', '36', '37', '38', '39','40',
                           '41', '42', '43', '44', '45'
                           ])
    hitting_times1 = mc.hitting_times(['1', '2', '3', '4', '5', '6', '7', '8', '9'])  # Target state C
    mean_hitting_times1 = np.mean(hitting_times1)
    #var_hitting_times1 = np.var(hitting_times1)

    hitting_times2 = mc.hitting_times(['10', '11', '12', '13', '14', '15', '16', '17', '18'])  # Target state C
    mean_hitting_times2 = np.mean(hitting_times2)
    #var_hitting_times2 = np.var(hitting_times2)

    hitting_times3 = mc.hitting_times(['19','20','21', '22', '23', '24', '25', '26', '27'])  # Target state C
    mean_hitting_times3 = np.mean(hitting_times3)
    #var_hitting_times3 = np.var(hitting_times3)

    hitting_times4 = mc.hitting_times(['28', '29','30', '31', '32', '33', '34', '35', '36'])  # Target state C
    mean_hitting_times4 = np.mean(hitting_times4)
    #var_hitting_times4 = np.var(hitting_times4)

    hitting_times5 = mc.hitting_times(['37', '38', '39','40','41', '42', '43', '44', '45'])  # Target state C
    mean_hitting_times5 = np.mean(hitting_times5)
    #var_hitting_times5 = np.var(hitting_times5)

    return [mean_hitting_times1, mean_hitting_times2, mean_hitting_times3, mean_hitting_times4, mean_hitting_times5]


def compute_ht2(data_array, log_edges):
    max_val = np.max(data_array)
    if max_val == 0:
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    log_edges_norm = (log_edges.max()-log_edges[::-1]) / (log_edges.max() - log_edges.min())
    scaled_edges = log_edges_norm * max_val

    if len(np.unique(scaled_edges)) != len(scaled_edges):
        return [0.0, 0.0, 0.0, 0.0, 0.0]

    num_states = len(scaled_edges) - 1
    transition_matrix = np.zeros((num_states, num_states))

    if len(data_array.shape) == 1:
        series = data_array
        states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
        states = pd.Series(states).fillna(0).astype(int)
        for i in range(len(states) - 1):
            current = states[i]
            next_ = states[i + 1]
            transition_matrix[current, next_] += 1
    else:
        for col in range(data_array.shape[1]):
            series = data_array[:, col]
            states = pd.cut(series, bins=scaled_edges, labels=False, include_lowest=True)
            states = pd.Series(states).fillna(0).astype(int)
            for i in range(len(states) - 1):
                current = states[i]
                next_ = states[i + 1]
                transition_matrix[current, next_] += 1

    row_sums = transition_matrix.sum(axis=1, keepdims=True)
    transition_probs = np.nan_to_num(transition_matrix / row_sums)
    mat = np.array(transition_probs)
    for i in range(mat.shape[0]):
        if mat[i].sum() == 0:
            mat[i] = np.full(mat.shape[1], 1.0 / mat.shape[1])

    mc = MarkovChain(mat, ['1', '2', '3', '4', '5', '6', '7', '8', '9','10', '11', '12', '13', '14', '15', '16', '17', '18', '19','20',
                           '21', '22', '23', '24', '25', '26', '27', '28', '29','30', '31', '32', '33', '34', '35', '36', '37', '38', '39','40',
                           '41', '42', '43', '44', '45'
                           ])
    hitting_times1 = mc.hitting_times(['1', '2', '3', '4', '5', '6', '7', '8', '9'])  # Target state C
    mean_hitting_times1 = np.mean(hitting_times1)
    #var_hitting_times1 = np.var(hitting_times1)

    hitting_times2 = mc.hitting_times(['10', '11', '12', '13', '14', '15', '16', '17', '18'])  # Target state C
    mean_hitting_times2 = np.mean(hitting_times2)
    #var_hitting_times2 = np.var(hitting_times2)

    hitting_times3 = mc.hitting_times(['19','20','21', '22', '23', '24', '25', '26', '27'])  # Target state C
    mean_hitting_times3 = np.mean(hitting_times3)
    #var_hitting_times3 = np.var(hitting_times3)

    hitting_times4 = mc.hitting_times(['28', '29','30', '31', '32', '33', '34', '35', '36'])  # Target state C
    mean_hitting_times4 = np.mean(hitting_times4)
    #var_hitting_times4 = np.var(hitting_times4)

    hitting_times5 = mc.hitting_times(['37', '38', '39','40','41', '42', '43', '44', '45'])  # Target state C
    mean_hitting_times5 = np.mean(hitting_times5)
    #var_hitting_times5 = np.var(hitting_times5)

    return [mean_hitting_times1, mean_hitting_times2, mean_hitting_times3, mean_hitting_times4, mean_hitting_times5]





##according to the chosen value of rho/theta, compare the distribution of psi and B1 along the sequence (100MB, window size of 1MB)
psi = []
psi2 = []
psi3 = []
B1 = []
B2 = []
#C = []
#Cv = []
Tpi = []
Wtheta = np.loadtxt(f"segSites_{sim_number}.txt", dtype=int).tolist()
n = 30  # Number of haploids
a_n = sum(1.0 / i for i in range(1, n))  # Harmonic number
num_cluster = 4

# Convert segregating sites to Watterson's theta by dividing by a_n
Wtheta = [s / a_n for s in Wtheta]

for i in range(len(distance_matrices)):
    psi.append(Psi(i, distance_matrices))
    psi2.append(Psi2(i, distance_matrices))
    psi3.append(Psi3(i, distance_matrices))
    B1.append(b1(i, distance_matrices))
    B2.append(b2(i, distance_matrices))
    #C.append(cluster(i, distance_matrices, num_cluster))
    #Cv.append(clusterv(i, distance_matrices, num_cluster))
    Tpi.append(tpi(i, distance_matrices, 30))

TajimasD = []
for i in range(len(distance_matrices)):
    D = tajimas_D(Tpi[i], int(Wtheta[i] * a_n), n)
    TajimasD.append(D)

mean_TajimasD = np.mean(TajimasD)
var_TajimasD = np.var(TajimasD)

##calculate the mean and variance of psi, B1, tajima's pi, and watterson's theta

mean_psi = np.mean(psi) ##
var_psi = np.var(psi) ##

mean_psi2 = np.mean(psi2) ##
var_psi2 = np.var(psi2) ##

mean_psi3 = np.mean(psi3) ##
var_psi3 = np.var(psi3) ##

mean_B1 = np.mean(B1) ##
var_B1 = np.var(B1) ##

mean_B2 = np.mean(B2) ##
var_B2 = np.var(B2) ##

#mean_C = np.mean(C) ##
#var_C = np.var(C) ##

#mean_Cv = np.mean(Cv) ##
#var_Cv = np.var(Cv) ##

mean_Tpi = np.mean(Tpi)
var_Tpi = np.var(Tpi)

mean_Wtheta = np.mean(Wtheta)
var_Wtheta = np.var(Wtheta)

delta1 = np.loadtxt(f"Delta1_{sim_number}.txt")
delta1_mean = np.mean(delta1, axis=1) ##along the window
delta1_mean2 = np.mean(delta1, axis=0) ##along the pairwise individuals
Delta1 = np.mean(delta1)
var_Delta1 = np.var(delta1_mean)
var2_Delta1 = np.var(delta1_mean2)

delta2 = np.loadtxt(f"Delta2_{sim_number}.txt")
delta2_mean = np.mean(delta2, axis=1)
delta2_mean2 = np.mean(delta2, axis=0)
Delta2 = np.mean(delta2)
var_Delta2 = np.var(delta2_mean)
var2_Delta2 = np.var(delta2_mean2)

delta3 = np.loadtxt(f"Delta3_{sim_number}.txt")
delta3_mean = np.mean(delta3, axis=1)
delta3_mean2 = np.mean(delta3, axis=0)
Delta3 = np.mean(delta3)
var_Delta3 = np.var(delta3_mean)
var2_Delta3 = np.var(delta3_mean2)

delta4 = np.loadtxt(f"Delta4_{sim_number}.txt",)
delta4_mean = np.mean(delta4, axis=1)
delta4_mean2 = np.mean(delta4, axis=0)
Delta4 = np.mean(delta4)
var_Delta4 = np.var(delta4_mean)
var2_Delta4 = np.var(delta4_mean2)

delta5 = np.loadtxt(f"Delta5_{sim_number}.txt")
delta5_mean = np.mean(delta5, axis=1)
delta5_mean2 = np.mean(delta5, axis=0)
Delta5 = np.mean(delta5)
var_Delta5 = np.var(delta5_mean)
var2_Delta5 = np.var(delta5_mean2)

delta6 = np.loadtxt(f"Delta6_{sim_number}.txt")
delta6_mean = np.mean(delta6, axis=1)
delta6_mean2 = np.mean(delta6, axis=0)
Delta6 = np.mean(delta6)
var_Delta6 = np.var(delta6_mean)
var2_Delta6 = np.var(delta6_mean2)

delta7 = np.loadtxt(f"Delta7_{sim_number}.txt")
delta7_mean = np.mean(delta7, axis=1)
delta7_mean2 = np.mean(delta7, axis=0)
Delta7 = np.mean(delta7)
var_Delta7 = np.var(delta7_mean)
var2_Delta7 = np.var(delta7_mean2)

delta8 = np.loadtxt(f"Delta8_{sim_number}.txt")
delta8_mean = np.mean(delta8, axis=1)
delta8_mean2 = np.mean(delta8, axis=0)
Delta8 = np.mean(delta8)
var_Delta8 = np.var(delta8_mean)
var2_Delta8 = np.var(delta8_mean2)

delta9 = np.loadtxt(f"Delta9_{sim_number}.txt")
delta9_mean = np.mean(delta9, axis=1)
delta9_mean2 = np.mean(delta9, axis=0)
Delta9 = np.mean(delta9)
var_Delta9 = np.var(delta9_mean)
var2_Delta9 = np.var(delta9_mean2)

re = np.loadtxt(f"r_{sim_number}.txt")
re_mean = np.mean(re, axis=1)
re_mean2 = np.mean(re, axis=0)
r_EMIBD = np.mean(re)
var_r_EMIBD = np.var(re_mean)
var2_r_EMIBD = np.var(re_mean2)

f_Mean = np.loadtxt(f"F_Mean_{sim_number}.txt")
f_Mean_mean = np.mean(f_Mean, axis=1)
f_Mean_mean2 = np.mean(f_Mean, axis=0) #along the individuals
F_Mean = np.mean(f_Mean)
var_F_Mean = np.var(f_Mean_mean)
var2_F_Mean = np.var(f_Mean_mean2)

f_SD = np.loadtxt(f"F_SD_{sim_number}.txt")
f_SD_mean = np.mean(f_SD, axis=1)
f_SD_mean2 = np.mean(f_SD, axis=0)
F_SD = np.mean(f_SD)
var_F_SD = np.var(f_SD_mean)
var2_F_SD = np.var(f_SD_mean2)

###obtaining transition matrix
log_edges = np.array([
    0.00000000, 0.0025,0.005, 0.0075, 0.01, 0.015, 0.02, 0.02531781, 0.05129329, 0.07796154,
    0.10536052, 0.13353139, 0.16251893, 0.19237189, 0.22314355, 0.25489225, 0.28768207, 0.32158362,
    0.35667494, 0.39304259, 0.43078292, 0.47000363, 0.51082562, 0.55338524, 0.59783700, 0.64435702,
    0.69314718, 0.74444047, 0.79850770, 0.85566611, 0.91629073, 0.98082925, 1.04982212, 1.12393010,
    1.20397280, 1.29098418, 1.38629436, 1.49165488, 1.60943791, 1.74296931, 1.89711998, 2.07944154,
    2.30258509, 2.59026717, 2.99573227, 3.68887945
])

Wtheta = np.array(Wtheta)
#mat_Wtheta = compute_transition_difference(Wtheta, log_edges)[0]
#ur_Wtheta = compute_transition_difference(Wtheta, log_edges)[1]
#ht12_mean_Wtheta = compute_transition_difference(Wtheta, log_edges)[2]
#ht12_var_Wtheta = compute_transition_difference(Wtheta, log_edges)[3]
#ht34_mean_Wtheta = compute_transition_difference(Wtheta, log_edges)[4]
#ht34_var_Wtheta = compute_transition_difference(Wtheta, log_edges)[5]
#ht5_mean_Wtheta = compute_transition_difference(Wtheta, log_edges)[6]
#ht5_var_Wtheta = compute_transition_difference(Wtheta, log_edges)[7]
#ht67_mean_Wtheta = compute_transition_difference(Wtheta, log_edges)[8]
#ht67_var_Wtheta = compute_transition_difference(Wtheta, log_edges)[9]
#ht89_mean_Wtheta = compute_transition_difference(Wtheta, log_edges)[10]
#ht89_var_Wtheta = compute_transition_difference(Wtheta, log_edges)[11]

psi = np.array(psi)
#ll_psi = compute_transition_difference(psi, log_edges)[0]
#ur_psi = compute_transition_difference(psi, log_edges)[1]
#ht12_mean_psi = compute_transition_difference(psi, log_edges)[2]
#ht12_var_psi = compute_transition_difference(psi, log_edges)[3]
#ht34_mean_psi = compute_transition_difference(psi, log_edges)[4]
#ht34_var_psi = compute_transition_difference(psi, log_edges)[5]
#ht5_mean_psi = compute_transition_difference(psi, log_edges)[6]
#ht5_var_psi = compute_transition_difference(psi, log_edges)[7]
#ht67_mean_psi = compute_transition_difference(psi, log_edges)[8]
#ht67_var_psi = compute_transition_difference(psi, log_edges)[9]
#ht89_mean_psi = compute_transition_difference(psi, log_edges)[10]
#ht89_var_psi = compute_transition_difference(psi, log_edges)[11]

psi2 = np.array(psi2)
#ll_psi2 = compute_transition_difference(psi2, log_edges)[0]
#ur_psi2 = compute_transition_difference(psi2, log_edges)[1]
#ht12_mean_psi2 = compute_ht(psi2, log_edges)[0]
#ht12_var_psi2 = compute_ht(psi2, log_edges)[1]
#ht34_mean_psi2 = compute_ht(psi2, log_edges)[2]
#ht34_var_psi2 = compute_ht(psi2, log_edges)[3]
#ht5_mean_psi2 = compute_ht(psi2, log_edges)[4]
#ht5_var_psi2 = compute_ht(psi2, log_edges)[5]
#ht67_mean_psi2 = compute_ht(psi2, log_edges)[6]
#ht67_var_psi2 = compute_ht(psi2, log_edges)[7]
#ht89_mean_psi2 = compute_transition_difference(psi2, log_edges)[10]
#ht89_var_psi2 = compute_transition_difference(psi2, log_edges)[11]

mat11 = []
mat12 = []
mat13 = []
mat14 = []
mat15 = []
mat21 = []
mat22 = []
mat23 = []
mat24 = []
mat25 = []
mat31 = []
mat32 = []
mat33 = []
mat34 = []
mat35 = []
mat41 = []
mat42 = []
mat43 = []
mat44 = []
mat45 = []
mat51 = []
mat52 = []
mat53 = []
mat54 = []
mat55 = []
var_lastcolum = []
var_lastrow = []
eigen_var = []
zero_count = []
ht1 = []
ht2 = []
ht3 = []
ht4 = []
ht5 = []

value2 = [Wtheta, psi, psi2]
value = [delta1, delta2, delta3, delta4, delta5, delta6, delta7, delta8, delta9, re, f_Mean, f_SD]

for i in value2:
    dif = compute_transition_difference(i, log_edges)
    dif2 = compute_ht(i, log_edges)
    mat11.append(dif[0])
    mat12.append(dif[1])
    mat13.append(dif[2])
    mat14.append(dif[3])
    mat15.append(dif[4])
    mat21.append(dif[5])
    mat22.append(dif[6])
    mat23.append(dif[7])
    mat24.append(dif[8])
    mat25.append(dif[9])
    mat31.append(dif[10])
    mat32.append(dif[11])
    mat33.append(dif[12])
    mat34.append(dif[13])
    mat35.append(dif[14])
    mat41.append(dif[15])
    mat42.append(dif[16])
    mat43.append(dif[17])
    mat44.append(dif[18])
    mat45.append(dif[19])
    mat51.append(dif[20])
    mat52.append(dif[21])
    mat53.append(dif[22])
    mat54.append(dif[23])
    mat55.append(dif[24])
    var_lastcolum.append(dif[25])
    var_lastrow.append(dif[26])
    eigen_var.append(dif[27])
    zero_count.append(dif[28])

    ##hitting times
    ht1.append(dif2[0])
    #var_ht1.append(dif2[1])
    ht2.append(dif2[1])
    #var_ht2.append(dif2[3])
    ht3.append(dif2[2])
    #var_ht3.append(dif2[5])
    ht4.append(dif2[3])
    #var_ht4.append(dif2[7])
    ht5.append(dif2[4])
    #var_ht5.append(dif2[9])

for i in value:
    dif = compute_transition_difference(i, log_edges)
    dif2 = compute_ht(i, log_edges)
    mat11.append(dif[0])
    mat12.append(dif[1])
    mat13.append(dif[2])
    mat14.append(dif[3])
    mat15.append(dif[4])
    mat21.append(dif[5])
    mat22.append(dif[6])
    mat23.append(dif[7])
    mat24.append(dif[8])
    mat25.append(dif[9])
    mat31.append(dif[10])
    mat32.append(dif[11])
    mat33.append(dif[12])
    mat34.append(dif[13])
    mat35.append(dif[14])
    mat41.append(dif[15])
    mat42.append(dif[16])
    mat43.append(dif[17])
    mat44.append(dif[18])
    mat45.append(dif[19])
    mat51.append(dif[20])
    mat52.append(dif[21])
    mat53.append(dif[22])
    mat54.append(dif[23])
    mat55.append(dif[24])
    var_lastcolum.append(dif[25])
    var_lastrow.append(dif[26])
    eigen_var.append(dif[27])
    zero_count.append(dif[28])

    ##hitting times
    ht1.append(dif2[0])
    #var_ht1.append(dif2[1])
    ht2.append(dif2[1])
    #var_ht2.append(dif2[3])
    ht3.append(dif2[2])
    #var_ht3.append(dif2[5])
    ht4.append(dif2[3])
    #var_ht4.append(dif2[7])
    ht5.append(dif2[4])
    #var_ht5.append(dif2[9])






# Write the full row to the CSV
with open("r1_reference_table7.csv", "a", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        sim_number,
        alpha,
        popsize,
        mean_psi, var_psi, 
        mean_psi2, var_psi2, 
        mean_psi3, var_psi3,
        mean_B1, var_B1,
        mean_B2, var_B2,
        mean_Tpi, var_Tpi,
        mean_Wtheta, var_Wtheta, 
        mean_TajimasD, var_TajimasD,
        Delta1, Delta2, Delta3, Delta4, Delta5, Delta6, Delta7, Delta8, Delta9, r_EMIBD, F_Mean, F_SD,
        var_Delta1, var_Delta2, var_Delta3, var_Delta4, var_Delta5, var_Delta6, var_Delta7, var_Delta8, var_Delta9, var_r_EMIBD, var_F_Mean, var_F_SD,
        var2_Delta1, var2_Delta2, var2_Delta3, var2_Delta4, var2_Delta5, var2_Delta6, var2_Delta7, var2_Delta8, var2_Delta9, var2_r_EMIBD, var2_F_Mean, var2_F_SD
    ] + mat11 + mat12 + mat13 + mat14 + mat15 + mat21 + mat22 + mat23 + mat24 + mat25 + mat31 + mat32 + mat33 + mat34 + mat35 + mat41 + mat42 + mat43 + mat44 + mat45 + mat51 + mat52 + mat53 + mat54 + mat55 + var_lastcolum + var_lastrow + eigen_var + zero_count + ht1 + ht2 + ht3 + ht4 + ht5)