import tskit
import msprime
from IPython.display import SVG
import numpy as np
import time
import matplotlib.pyplot as plt
from collections import Counter
import sys
import csv
import subprocess
import scipy
import sgtl
import sgtl.graph
import sgtl.spectrum
import sgtl.random
import sgtl.clustering
import scipy.sparse as sp
import pandas as pd
from pydtmc import MarkovChain
import itertools as itertools

#genealogy simuation
seq_length = 5_000_000 ##1MB
Ttime = [0, 100, 200, 5e2, 1000, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5]
transition_time = int(np.random.choice(Ttime))

#print(f"Simulating with transition_time = {transition_time}")
alpha = 1.3
#print(f"Simulating with alpha = {alpha}")
r = 1e-8
popsize = 5e4
#print(f"Simulating with popsize = {popsize}")
#to simulate the sequence data, we need to define the proportion of sites of each marker type in the sequence
#need to make sure if the logic is correct

ts = msprime.sim_ancestry(
    15,
    recombination_rate=r,
    sequence_length=seq_length,
    population_size=popsize,
    model=[
        msprime.StandardCoalescent(duration = transition_time),
        msprime.BetaCoalescent(alpha=alpha)
    ]
)

samples = list(ts.samples())
TMRCA = []
for a,b in itertools.combinations(samples,2):
    tmrca = []
    for tree in ts.trees():
        mrca_node = tree.mrca(a, b)
        t = tree.time(mrca_node)
        tmrca.append(t)
    TMRCA.append(tmrca)
tmrca_matrix = np.array(TMRCA)  # ex Shape: (435, 411)
tmrca_matrix = tmrca_matrix.T  # ex Shape: (411, 435)

def compute_transition_difference(data_array, log_edges):
    max_val = np.max(data_array)
    if max_val == 0:
        return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0,0.0,0.0]

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

log_edges = np.array([
    0.00000000, 0.0025,0.005, 0.0075, 0.01, 0.015, 0.02, 0.02531781, 0.05129329, 0.07796154,
    0.10536052, 0.13353139, 0.16251893, 0.19237189, 0.22314355, 0.25489225, 0.28768207, 0.32158362,
    0.35667494, 0.39304259, 0.43078292, 0.47000363, 0.51082562, 0.55338524, 0.59783700, 0.64435702,
    0.69314718, 0.74444047, 0.79850770, 0.85566611, 0.91629073, 0.98082925, 1.04982212, 1.12393010,
    1.20397280, 1.29098418, 1.38629436, 1.49165488, 1.60943791, 1.74296931, 1.89711998, 2.07944154,
    2.30258509, 2.59026717, 2.99573227, 3.68887945
])

compartments = compute_transition_difference(tmrca_matrix, log_edges)


# Write the full row to the CSV
with open("r1_reference_table7.csv", "a", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        transition_time,
        alpha,
        popsize,
    ] + compartments)