import matplotlib.pyplot as plt # for plotting

import numpy as np
import pulp
from pulp import *

import scipy # various algorithms
import scipy.optimize
import scipy.interpolate
import scipy.signal
import random
import pandas as pd
import scipy as sp

import time

solverdir = "C:\Users\Nicholas\Anaconda2\Lib\site-packages\pulp\libCbcSolver.so"
solver = pulp.COIN_CMD(path=solverdir)
# Load data

df = pd.read_csv('C:/Users/Nicholas/OneDrive/Berkeley Biophysics/projects/all_the_stripes/data/all_traces.csv', header=None)
#print(pd.DataFrame.head(df))

rows = [i for i in range(0, len(df))]
rows = [21]

fig = plt.figure(figsize=(8, 4 * len(rows)))
subplot = 0

# determined empirically

pulse_duration = 16
pulse_intensity = 25000

for d in rows:

    # data = np.nan_to_num(np.array(df.loc[d][1:]))  # first column is AP position - don't use for modeling
    data = np.array(df.iloc[d, :])
    datalen = len(data)

    # just for convenience later
    data_indices = [i for i in range(0, len(data))]

    # create basis vectors corresponding to pulses

    basis_vectors = []
    for x in range(0, datalen):
        basis_vector = datalen * [0]
        for i in range(x, min(x + pulse_duration, datalen)):
            basis_vector[i] = pulse_intensity
        basis_vectors.append(np.array(basis_vector))

    # for convenience
    basis_indices = [i for i in range(len(basis_vectors))]

    # create variable for PuLP

    # basis loadings is 0 or 1 corresponding to if basis vector is used
    basis_loadings = [LpVariable("Basis vector " + str(b), 0, 1, cat="Integer") for b in basis_indices]

    # variable to capture to error in estimates
    err = [LpVariable("err " + str(i), 0, None, cat="Continuous") for i in data_indices]

    # create a PuLP problem - going to minimize objective function
    prob = LpProblem("Fit MS2 Pulses", LpMinimize)

    # objective function is sum of errors
    prob += sum([err[i] for i in data_indices])

    # constraint is that modeled data less than error and greater than -error
    for i in data_indices:
        prob += data[i] - sum([basis_loadings[j] * basis_vectors[j][i] for j in basis_indices]) <= err[i]
        prob += data[i] - sum([basis_loadings[j] * basis_vectors[j][i] for j in basis_indices]) >= -err[i]

    # The problem data is written to an .lp file
    prob.writeLP("BasisTest.lp")

    # Solve problem
    #prob.solve(pulp.GLPK())
    prob.solve()

    # Calculate fit data

    fitdata = datalen * [0]
    fitdata = np.array(fitdata)
    basis_loading_vector = []

    for b in basis_indices:
        if basis_loadings[b].varValue == 1.0:
            fitdata += basis_vectors[b]
        basis_loading_vector.append(basis_loadings[b].varValue)
    # Plot

    subplot += 1
    ax = plt.subplot(len(rows), 1, subplot)
    ax.plot(data, alpha=0.4)
    ax.plot(fitdata, color='r')
    ax.set_title("Trace " + str(d) + " with pulse intensity " + str(pulse_intensity) + " and pulse duration " + str(
        pulse_duration))
    ax2 = ax.twinx()
    ax2.bar(np.arange(0, len(basis_loading_vector)), basis_loading_vector, color='g', edgecolor="none", alpha=0.05)

plt.tight_layout()
plt.show()