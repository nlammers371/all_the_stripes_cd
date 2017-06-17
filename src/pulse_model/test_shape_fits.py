import scipy # various algorithms
import scipy.optimize
import scipy.interpolate
import scipy.signal
import scipy.linalg
import random
import pandas as pd
import scipy as sp
import os
import time
import numpy as np
from matplotlib import pyplot as plt
#Import Trace Generator
from generate_traces import generate_traces_gill, generate_traces_prob_mat
from functions import shape_mode_decomposition, shape_mode_fits
#Import Model(s)
from models import get_pulse_sequence
from itertools import chain
from models_alternative import get_pulse_sequence_continuous, get_pulse_sequence_multi
#-----------------------------------------------Test Parameters--------------------------------------------------------#
#Num Traces
n_traces = 4
#Trace length (in time steps)
tr_length = 160
#Set Memory (in time steps)
w = 16
#Set noise scale (relative to max possible fluo)
noise = .05
#Set rise time for MS2
alpha = 0
#Transition rate matrix (rates in timestep units)
R = np.array([[-.01, .007, .00],
              [.01, -.015, .02],
              [.000, .008, -.02]]) * 10.2

#R = np.array([[-.008, .015], [.008, -.015]]) * 10.2
#Transition Prob matrix
A = scipy.linalg.expm(R, q=None)

#Initiation rates (fluo units)
#v = np.array([0, 25])
v = np.array([0, 25, 50])
#Trace type
ideal_traces = False
#Set writepath for results
outpath = '../../results/shape_fits'
#Set project name (creates subfolder)
project_name = 'comparison_test'
if not os.path.isdir(os.path.join(outpath,project_name)):
    os.makedirs(os.path.join(outpath,project_name))

#----------------------------------------Continuous Model Params-------------------------------------------------------#
cont=1
#Minimum unit of variation in traces (
r_unit = v[1]
#if 0, discrete penalty in objective is ignored
discrete=1
#Set scale of diff penalty
inertia = 1.5
#------------------------------------------Discrete Model Params-------------------------------------------------------#
unit_pulse = 0
pulse_size = v[1]
n_states = 3
#-----------------------------------------------Shape Fit Params-------------------------------------------------------#
n_modes = 10 #Must be <= w
#Num Traces used to generate shape eigenvectors
n_gen_traces = 500
#---------------------------------------Generate Shape Vectors---------------------------------------------------------#
eigenvalues, eigenvectors, fragment_mean = shape_mode_decomposition(memory=w,
                                                                    length=tr_length,
                                                                    input_size=w*max(v)+1,
                                                                    batch_size=n_gen_traces,
                                                                    r_mat=R,
                                                                    v=v,
                                                                    noise_scale=.025,
                                                                    alpha=1.0)
#--------------------------------------------Generate Traces-----------------------------------------------------------#

if ideal_traces:
    fluo_seq, int_labels= generate_traces_prob_mat(A=A,
                                                      v=v,
                                                      memory=w,
                                                      length=tr_length,
                                                      input_size=w*max(v)+1,
                                                      batch_size=n_traces,
                                                      alpha=alpha,
                                                      noise_scale=noise)

else:
    fluo_mat, pulse_mat, seq_lengths, int_labels, fluo_seq, true_seq = generate_traces_gill(memory=w,
                                                                                                length=tr_length,
                                                                                                input_size=w*max(v)+1,
                                                                                                batch_size=n_traces,
                                                                                                r_mat=R,
                                                                                                v=v,
                                                                                                noise_scale = noise,
                                                                                                alpha=alpha)
#--------------------------------------------Generate Shape Fits-------------------------------------------------------#
traces = np.array(fluo_seq)
full_mat, shape_fits = shape_mode_fits(np.real(eigenvectors), fragment_mean, np.matrix.transpose(traces), n_modes)

#---------------------------------------------Run Fitting Routine------------------------------------------------------#
if unit_pulse:
    f_solutions_unit = []
    p_solutions_unit = []
if cont:
    f_solutions_cont = []
    p_solutions_cont = []

fluo_true_states = []
fluo_seq_out = []
int_labels_out = []

for i in xrange(2*len(fluo_seq)):
    if i % 2 == 0:
        fluo_seq_out.append(fluo_seq[i/2])
        fluo_true_states.append(true_seq[i/2])
    else:
        fluo_seq_out.append(shape_fits[i/2])
        fluo_true_states.append(true_seq[i / 2])
    int_labels_out.append(int_labels[i / 2])
#add shape fits to f vec list

#print(int_labels)
for tr in xrange(len(fluo_seq_out)):

    if unit_pulse:
        print("Conducting Unit Pulse Fit " + str(tr) + "...")
        init_time = time.time()
        basis_loading_vector_unit, fitdata_unit = get_pulse_sequence_multi(fluo_seq_out[tr], num_states=n_states, smooth=True, pulse_intensity=pulse_size, pulse_duration=w, smoothp=2*v[1])
        basis_loading_vector_unit = pulse_size * basis_loading_vector_unit
        p_solutions_unit.append(basis_loading_vector_unit)
        f_solutions_unit.append(fitdata_unit)

        fig_fluo = plt.figure(figsize=(12, 4))
        ax = plt.subplot(1, 1, 1)
        ax.plot(fluo_seq_out[tr], c='b', alpha=0.4, label='Actual')
        ax.plot(f_solutions_unit[tr], c='r', label='Predicted')
        # ax.plot(newdata,c='g')
        # plt.show()
        #plt.legend()
        fig_fluo.savefig(os.path.join(outpath, project_name, 'f' "_" + str(tr/2) + '_' + str(1 + tr % 2) + '_unit_trajectory' +  ".png"))
        plt.close()

        fig_prom = plt.figure(figsize=(12, 4))
        ax = plt.subplot(1, 1, 1)
        ax.plot(int_labels_out[tr], c='b', alpha=0.4, label='Actual')
        ax.plot(p_solutions_unit[tr], c='r', label='Predicted')
        # ax.plot(newdata,c='g')
        # plt.show()
        #plt.legend()
        fig_prom.savefig(os.path.join(outpath, project_name, 'p' + "_" + str(tr/2) + '_' + str(1 + tr % 2) + 'unit_trajectory'  ".png"))
        plt.close()
        print("Runtime (unit): " + str(time.time() - init_time))

    if cont:
        print("Conducting Continuous Fit " + str(tr) + "...")
        init_time = time.time()
        init_diffs = np.diff(fluo_seq_out[tr])
        m_diff = np.max(np.abs(init_diffs))
        m_diff_diff = np.max(np.abs(np.diff(init_diffs, 2)))

        basis_loading_vector_cont, fitdata_cont = get_pulse_sequence_continuous(data=fluo_seq_out[tr], diff_init=init_diffs, max_diff = m_diff_diff, u_bound=m_diff, inertia = inertia, discrete=0, r_unit=r_unit, pulse_duration=w)
        p_solutions_cont.append(basis_loading_vector_cont)
        f_solutions_cont.append(fitdata_cont)

        fig_fluo = plt.figure(figsize=(12, 4))
        ax = plt.subplot(1, 1, 1)
        ax.plot(fluo_seq_out[tr], c='b', alpha=0.4, label='Actual')
        ax.plot(f_solutions_cont[tr], c='r', label='Predicted')
        # ax.plot(newdata,c='g')
        # plt.show()
        #plt.legend()
        fig_fluo.savefig(os.path.join(outpath, project_name, 'f' + "_" + str(tr/2) + '_' + str(1 + tr % 2) + "_cont_trajectory.png"))
        plt.close()

        fig_prom = plt.figure(figsize=(12, 4))
        ax = plt.subplot(1, 1, 1)
        ax.plot(int_labels_out[tr], c='b', alpha=0.4, label='Actual')
        ax.plot(p_solutions_cont[tr], c='r', label='Predicted')
        # ax.plot(newdata,c='g')
        # plt.show()
        #plt.legend()
        fig_prom.savefig(os.path.join(outpath, project_name, 'p' + "_" + str(tr/2) + '_' + str(1 + tr % 2) + '_cont_trajectory.png'))
        plt.close()
        print("Runtime (cont): " + str(time.time() - init_time))

    fig_fluo = plt.figure(figsize=(12, 4))
    ax = plt.subplot(1, 1, 1)
    ax.plot(fluo_seq_out[tr], c='b', alpha=0.4, label='Appx')
    ax.plot(fluo_true_states[tr], c='r', label='Actual')
    #ax.plot(fluo_true_states[tr], c='r', label='Predicted')
    # ax.plot(newdata,c='g')
    # plt.show()
    # plt.legend()
    fig_fluo.savefig(os.path.join(outpath, project_name, 'f' + "_" + str(tr / 2) + '_' + str(1 + tr % 2) + "_comparison.png"))
    plt.close()


