import numpy as np
from pulp import *

# this will fit a pulse sequence from data
# if smooth is False it will just fit sequence
# if smooth is True it will apply smoothing, but runs slowly

def get_pulse_sequence(data, smooth=False, pulse_intensity=25000, pulse_duration=16, smoothp=100000):
    data = np.nan_to_num(data)
    datalen = len(data)
    data_indices = [i for i in range(0, datalen)]

    basis_vectors = []
    for x in range(0, datalen):
        basis_vector = datalen * [0]
        for i in range(x, min(x + pulse_duration, datalen)):
            basis_vector[i] = pulse_intensity
        basis_vectors.append(np.array(basis_vector))

    basis_indices = [i for i in range(len(basis_vectors))]

    basis_loadings = [LpVariable("Basis vector " + str(b), 0, 1, cat="Integer") for b in basis_indices]
    err = [LpVariable("err " + str(i), 0, None, cat="Continuous") for i in data_indices]

    if smooth == True:
        loadshifts = []
        for b in basis_indices[0:-1]:
            loadshifts.append(LpVariable("Loadshift " + str(b), 0, 1, cat="Continuous"))

    prob = LpProblem("Fit MS2 Pulses", LpMinimize)

    if smooth == True:
        prob += sum([err[i] for i in data_indices]) + smoothp * sum([loadshifts[b] for b in basis_indices[0:-1]])
    else:
        prob += sum([err[i] for i in data_indices])

    for i in data_indices:
        prob += data[i] - sum([basis_loadings[j] * basis_vectors[j][i] for j in basis_indices]) <= err[i]
        prob += data[i] - sum([basis_loadings[j] * basis_vectors[j][i] for j in basis_indices]) >= -err[i]

    if smooth == True:
        for b in basis_indices[0:-1]:
            prob += basis_loadings[b + 1] - basis_loadings[b] <= loadshifts[b]
            prob += basis_loadings[b + 1] - basis_loadings[b] >= -loadshifts[b]

    prob.writeLP("MS2Fit.lp")
    prob.solve()

    fitdata = datalen * [0]
    fitdata = np.array(fitdata)

    basis_loading_vector = np.zeros(len(data))

    for b in basis_indices:
        if basis_loadings[b].varValue == 1.0:
            fitdata += basis_vectors[b]
            basis_loading_vector[b] = 1

    return (basis_loading_vector, fitdata)



# to work around slowness of smoothing
# this function fixes all of the pulses except in a specified window
# refits in that window only using smoothing penalty

def refine_pulse_sequence(data, pulse_sequence, win_start, win_end, pulse_intensity=25000, pulse_duration=16,
                          smoothp=100000):
    data = np.nan_to_num(data)
    datalen = len(data)
    data_indices = [i for i in range(0, datalen)]

    basis_vectors = []
    for x in range(0, datalen):
        basis_vector = datalen * [0]
        for i in range(x, min(x + pulse_duration, datalen)):
            basis_vector[i] = pulse_intensity
        basis_vectors.append(np.array(basis_vector))

    basis_indices = [i for i in range(len(basis_vectors))]

    basis_loadings = []

    for b in basis_indices:
        if b in range(win_start, win_end):
            basis_loadings.append(LpVariable("Basis vector " + str(b), 0, 1, cat="Integer"))
        else:
            basis_loadings.append(
                LpVariable("Basis vector " + str(b), pulse_sequence[b], pulse_sequence[b], cat="Integer"))

    err = [LpVariable("err " + str(i), 0, None, cat="Continuous") for i in data_indices]

    loadshifts = (datalen - 1) * [0]
    for b in range(win_start, win_end):
        loadshifts[b] = LpVariable("Loadshift " + str(b), 0, 1, cat="Continuous")

    prob = LpProblem("Fit MS2 Pulses", LpMinimize)

    prob += sum([err[i] for i in data_indices]) + smoothp * sum([loadshifts[b] for b in range(win_start, win_end)])

    for i in data_indices:
        prob += data[i] - sum([basis_loadings[j] * basis_vectors[j][i] for j in basis_indices]) <= err[i]
        prob += data[i] - sum([basis_loadings[j] * basis_vectors[j][i] for j in basis_indices]) >= -err[i]

    for b in range(win_start, win_end - 1):
        prob += basis_loadings[b + 1] - basis_loadings[b] <= loadshifts[b]
        prob += basis_loadings[b + 1] - basis_loadings[b] >= -loadshifts[b]

    prob.writeLP("MS2Fit.lp")
    prob.solve()

    fitdata = datalen * [0]
    fitdata = np.array(fitdata)

    basis_loading_vector = np.zeros(len(data))

    for b in basis_indices:
        if basis_loadings[b].varValue == 1.0:
            fitdata += basis_vectors[b]
            basis_loading_vector[b] = 1

    return (np.array(basis_loading_vector), fitdata)

