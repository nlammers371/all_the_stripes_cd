from scipy.optimize import least_squares
import numpy as np
from pulp import *

#Continuouse fromulation of original optimization problem

def get_pulse_sequence_continuous(data, diff_init, max_diff, r_unit, inertia, u_bound, discrete=1.0, pulse_duration=16):
    data = np.nan_to_num(data)
    datalen = len(data)

    #Define array containing basis vectors
    basis_vectors = np.zeros((datalen,datalen), dtype='float')
    for x in range(0, datalen):
        basis_vectors[x,x:min(x+pulse_duration,datalen)] = 1.0

    def model(x, b):
        return np.matmul(x,b)
    #Define objective function
    def fun(x, b, f_data):
        return np.sum(np.square(f_data-model(x,b))) + \
                np.sum(np.square(inertia*(np.diff(x,2))))\
               + discrete*np.sum(np.square(4*r_unit*np.sin(np.pi*(x / r_unit))))

    #Define initialization vector for states
    x0 = np.zeros((datalen))
    x0[0] =0.0
    for i in xrange(1,datalen):
        x0[i] = max(0.0,diff_init[i-1])

    res = least_squares(fun, x0, bounds=(0, u_bound), args=(basis_vectors, data), verbose=1)

    basis_loading_vector = res.x

    fitdata = np.matmul(basis_loading_vector, basis_vectors)

    return (basis_loading_vector, fitdata)


def get_pulse_sequence_multi(data, num_states, smooth=False, pulse_intensity=25000, pulse_duration=16, smoothp=100000):
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

    basis_loadings = [LpVariable("Basis vector " + str(b), 0, num_states, cat="Integer") for b in basis_indices]
    err = [LpVariable("err " + str(i), 0, None, cat="Continuous") for i in data_indices]

    if smooth == True:
        loadshifts = []
        for b in basis_indices[0:-1]:
            loadshifts.append(LpVariable("Loadshift " + str(b), 0, num_states, cat="Continuous"))

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

    fitdata = np.zeros(datalen,dtype='int')

    basis_loading_vector = np.zeros(len(data))

    for b in basis_indices:
        fitdata += basis_vectors[b] * int(basis_loadings[b].varValue)
        basis_loading_vector[b] = basis_loadings[b].varValue

    return (basis_loading_vector, fitdata)
