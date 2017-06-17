import scipy # various algorithms
import scipy.optimize
import scipy.interpolate
import scipy.signal
import scipy.linalg
import random
import pandas as pd
import scipy as sp
from itertools import chain
import os
import time
import numpy as np
from matplotlib import pyplot as plt
#Import Trace Generator
from generate_traces import generate_traces_gill, generate_traces_prob_mat
#Import Model(s)
from models import get_pulse_sequence
from models_alternative import get_pulse_sequence_continuous, get_pulse_sequence_multi
from functions import on_off_streaks, get_pulse_size

#-------------------------------------------Global Test Parameters-----------------------------------------------------#
####IS THIS A SIMULATION????############
simulation = 1
#Num Traces
n_traces = 2
#Trace length (in time steps)
tr_length = 150
#Set Memory (in time steps)
w = 16
#Set rise time for MS2
alpha = w * 60 / 204
#Trace type
ideal_traces = False

#-------------------------------------Routine-Specific Test Parameters-------------------------------------------------#
#Set noise scale (relative to avg of v)
noise_list = [1.0, 1.0]

#Transition rate matrix (rates in timestep units)
#R_list = [np.array([[-.008, .007, .00],[.008, -.015, .07],[.000, .008, -.07]]) * 6.0,
#     np.array([[-.008, .007, .00], [.008, -.015, .04], [.000, .008, -.04]]) * 10.2]
R_list = [np.array([[-.008, .015],[.008, -.015]]) * 6.0,
     np.array([[-.008, .015], [.008, -.015]]) * 10.2]
#Initiation rates (fluo units)
#v = np.array([0, 25])
#v_list = [np.array([0, 25, 50]),np.array([0, 25, 50])]
v_list = [np.array([0, 25]),np.array([0, 25])]

subfolder_names = ['fast_slow', 'fast_slow']
test_names = ['_2_state_slow' , '_2_state_fast']

##########Continuous Model Params#################
cont=0
#Minimum unit of variation in traces (
r_unit = 5.0
#if 0, discrete penalty in objective is ignored
discrete=1
#Set scale of diff penalty
inertia = 1.0
########Unit Model Params################
unit_pulse = 1
pulse_size = v_list[0][1]
n_on_states = 1

#Iterate Through Test Routines
for t in xrange(len(test_names)):
    R = R_list[t]
    v = v_list[t]
    noise = noise_list[t]
    #Transition Prob matrix
    A = scipy.linalg.expm(R, q=None)
    #Set writepath for results
    outpath = '../../results/unit_pulse_validation/' + subfolder_names[t]
    #Set project name (creates subfolder)
    subfolder_name = 'unit_pulse' + test_names[t]
    if not os.path.isdir(os.path.join(outpath,subfolder_name)):
        os.makedirs(os.path.join(outpath,subfolder_name))
    if not os.path.isdir(os.path.join(outpath,subfolder_name,'plots')):
        os.makedirs(os.path.join(outpath,subfolder_name,'plots'))
    print("Writing to: " +  os.path.join(outpath,subfolder_name))
    #--------------------------------------------Generate Traces-----------------------------------------------------------#
    if simulation:
        if ideal_traces:
            fluo_seq, int_labels = generate_traces_prob_mat(A=A,
                                                              v=v,
                                                              memory=w,
                                                              length=tr_length,
                                                              input_size=w*max(v)+1,
                                                              batch_size=n_traces,
                                                              alpha=alpha,
                                                              noise_scale=noise)

        else:
            fluo_mat, pulse_mat, seq_lengths, int_labels, fluo_seq, conv_labels = generate_traces_gill(memory=w,
                                                                                                        length=tr_length,
                                                                                                        input_size=w*max(v)+1,
                                                                                                        batch_size=n_traces,
                                                                                                        r_mat=R,
                                                                                                        v=v,
                                                                                                        noise_scale = noise,
                                                                                                        alpha=alpha)
    #else:

    #---------------------------------------------Run Fitting Routine------------------------------------------------------#
    if unit_pulse:
        f_solutions_unit = []
        p_solutions_unit = []
        off_lengths_unit = []
        on_lengths_unit = []
        pulse_sizes_unit = []
    if cont:
        f_solutions_cont = []
        p_solutions_cont = []
        off_lengths_cont = []
        on_lengths_cont = []
        pulse_sizes_cont = []

    off_lengths_true = []
    on_lengths_true = []
    pulse_sizes_true = []


    for tr in xrange(n_traces):
        if simulation:
            z_inf, nz_inf, change_points = on_off_streaks(np.array(int_labels[tr]))
            pulse_sizes_inf = get_pulse_size(np.array(int_labels[tr]), change_points)

            off_lengths_true.append(z_inf)
            on_lengths_true.append(nz_inf)
            pulse_sizes_true.append(pulse_sizes_inf)

        if unit_pulse:
            print("Conducting Unit Pulse Fit " + str(tr) + "...")
            init_time = time.time()
            basis_loading_vector_unit, fitdata_unit = get_pulse_sequence_multi(fluo_seq[tr], num_states=n_on_states, smooth=True, pulse_intensity=pulse_size, pulse_duration=w, smoothp=2*v[1])
            basis_loading_vector_unit = pulse_size * basis_loading_vector_unit

            z_inf, nz_inf, change_points = on_off_streaks(basis_loading_vector_unit)
            pulse_sizes_inf = get_pulse_size(basis_loading_vector_unit, change_points)

            off_lengths_unit.append(z_inf)
            on_lengths_unit.append(nz_inf)
            pulse_sizes_unit.append(pulse_sizes_inf)

            p_solutions_unit.append(basis_loading_vector_unit)
            f_solutions_unit.append(fitdata_unit)

            fig_fluo = plt.figure(figsize=(12, 4))
            ax = plt.subplot(1, 1, 1)
            if simulation:
                ax.plot(fluo_seq[tr], c='b', alpha=0.4, label='Actual')
            ax.plot(f_solutions_unit[tr], c='r', label='Predicted')

            #plt.legend()
            fig_fluo.savefig(os.path.join(outpath, subfolder_name, 'plots', 'f_trajectory_unit' + "_" + str(tr) + ".png"))
            plt.close()

            fig_prom = plt.figure(figsize=(12, 4))
            ax = plt.subplot(1, 1, 1)
            if simulation:
                ax.plot(int_labels[tr], c='b', alpha=0.4, label='Actual')
            ax.plot(p_solutions_unit[tr], c='r', label='Predicted')

            #plt.legend()
            fig_prom.savefig(os.path.join(outpath, subfolder_name, 'plots', 'p_trajectory_unit' + "_" + str(tr) + ".png"))
            plt.close()
            print("Runtime (unit): " + str(time.time() - init_time))



        if cont:
            print("Conducting Continuous Fit " + str(tr) + "...")
            init_time = time.time()
            init_diffs = np.diff(fluo_seq[tr])
            m_diff = np.max(np.abs(init_diffs))
            m_diff_diff = np.max(np.abs(np.diff(init_diffs, 2)))

            basis_loading_vector_cont, fitdata_cont = get_pulse_sequence_continuous(data=fluo_seq[tr], diff_init=init_diffs, max_diff = m_diff_diff, u_bound=m_diff, inertia = inertia, discrete=0, r_unit=r_unit, pulse_duration=w)

            z_inf, nz_inf, change_points = on_off_streaks(basis_loading_vector_cont)
            pulse_sizes_inf = get_pulse_size(basis_loading_vector_cont, change_points)

            off_lengths_cont.append(z_inf)
            on_lengths_cont.append(nz_inf)
            pulse_sizes_cont.append(pulse_sizes_inf)

            p_solutions_cont.append(basis_loading_vector_cont)
            f_solutions_cont.append(fitdata_cont)

            fig_fluo = plt.figure(figsize=(12, 4))
            ax = plt.subplot(1, 1, 1)
            if simulation:
                ax.plot(fluo_seq[tr], c='b', alpha=0.4, label='Actual')
            ax.plot(f_solutions_cont[tr], c='r', label='Predicted')

            #plt.legend()
            fig_fluo.savefig(os.path.join(outpath, subfolder_name, 'plots', 'f_trajectory_cont' + "_" + str(tr) + ".png"))
            plt.close()

            fig_prom = plt.figure(figsize=(12, 4))
            ax = plt.subplot(1, 1, 1)
            if simulation:
                ax.plot(int_labels[tr], c='b', alpha=0.4, label='Actual')
            ax.plot(p_solutions_cont[tr], c='r', label='Predicted')
            #plt.legend()
            fig_prom.savefig(os.path.join(outpath, subfolder_name, 'plots', 'p_trajectory_cont' + "_" + str(tr) + ".png"))
            plt.close()
            print("Runtime (cont): " + str(time.time() - init_time))

        if max(unit_pulse,cont) == 0 and simulation==1:
            print("Saving Sample Traces (No fits requested)")
            fig_fluo = plt.figure(figsize=(12, 4))
            ax = plt.subplot(1, 1, 1)
            ax.plot(fluo_seq[tr], c='b', alpha=0.4, label='Actual')
            # plt.legend()
            fig_fluo.savefig(os.path.join(outpath, subfolder_name, 'plots', 'f_trajectory_unit' + "_" + str(tr) + ".png"))
            plt.close()
            fig_prom = plt.figure(figsize=(12, 4))
            ax = plt.subplot(1, 1, 1)
            ax.plot(int_labels[tr], c='b', alpha=0.4, label='Actual')
            fig_prom.savefig(os.path.join(outpath, subfolder_name, 'plots', 'p_trajectory_cont' + "_" + str(tr) + ".png"))
    #---------------------------------------------------Make Histograms----------------------------------------------------#
    if max(unit_pulse,cont) > 0:
        #Off Times
        fig_off = plt.figure(figsize=(15,5))
        ax1 = fig_off.add_subplot(111)
        # the histogram of the data

        bins = np.linspace(0, 100, 100)
        if simulation:
            ax1.hist(np.array(list(chain(*off_lengths_true))), bins=bins, normed =1,edgecolor='black', alpha=0.5, label="Actual")
        if unit_pulse:
            ax1.hist(np.array(list(chain(*off_lengths_unit))), bins=bins, normed =1,edgecolor='black', alpha=0.5, label="Predicted (Unit)")
        if cont:
            ax1.hist(np.array(list(chain(*off_lengths_unit))), bins=bins, normed=1, edgecolor='black', alpha=0.5, label="Predicted (Continuous)")

        plt.xlabel('Duration (Time Steps)')
        plt.ylabel('Share')
        plt.title('Distribution of Off Times')
        plt.legend()
        plt.grid(True)
        fig_off.savefig(os.path.join(outpath, subfolder_name, 'off_distributions.png'))
        #plt.show()

        #On Times
        fig_on = plt.figure(figsize=(15,5))
        ax2 = fig_on.add_subplot(111)
        # the histogram of the data

        bins = np.linspace(0, 100, 100)
        if simulation:
            ax2.hist(np.array(list(chain(*on_lengths_true))), bins=bins, normed =1,edgecolor='black', alpha=0.5, label="Actual")
        if unit_pulse:
            ax2.hist(np.array(list(chain(*on_lengths_unit))), bins=bins, normed =1,edgecolor='black', alpha=0.5, label="Predicted (Unit)")
        if cont:
            ax2.hist(np.array(list(chain(*on_lengths_unit))), bins=bins, normed=1, edgecolor='black', alpha=0.5, label="Predicted (Continuous)")

        plt.xlabel('Duration (Time Steps)')
        plt.ylabel('Share')
        plt.title('Distribution of On Times')
        plt.legend()
        plt.grid(True)
        fig_on.savefig(os.path.join(outpath, subfolder_name, 'on_distributions.png'))
        #plt.show()

        #Pulse Sizes
        fig_size = plt.figure(figsize=(15,5))
        ax3 = fig_size.add_subplot(111)
        # the histogram of the data

        bins = np.linspace(0, 3000, 100)
        if simulation:
            ax3.hist(np.array(list(chain(*pulse_sizes_true))), bins=bins, normed =1,edgecolor='black', alpha=0.5, label="Actual")
        if unit_pulse:
            ax3.hist(np.array(list(chain(*pulse_sizes_unit))), bins=bins, normed=1,edgecolor='black', alpha=0.5, label="Predicted (Unit)")
        if cont:
            ax3.hist(np.array(list(chain(*pulse_sizes_unit))), bins=bins, normed=1, edgecolor='black', alpha=0.5, label="Predicted (Continuous)")

        plt.xlabel('Total Fluorescence (A.U.)')
        plt.ylabel('Share')
        plt.title('Distribution of Pulse Sizes')
        plt.legend()
        plt.grid(True)
        fig_size.savefig(os.path.join(outpath, subfolder_name, 'pulse_size_distributions.png'))
        #plt.show()

        #------------------------------------------------Make Scatter Plots----------------------------------------------------#

        #Pulse length visualization
        fig_pulse_scatter = plt.figure(figsize=(15,5))
        ax4 = fig_pulse_scatter.add_subplot(111)
        hand = []
        for i in xrange(min(50,len(pulse_sizes_true))):
            lb_axis = np.where(np.array(int_labels[i]) > 0)[0]
            if simulation==1:
                true_bursts = ax4.plot(lb_axis, np.linspace(3*i + 1,3*i + 1,len(lb_axis)),'o',c='green', alpha=.5, label = 'Actual')
                if i==0:
                    hand.append(true_bursts[0])
            if unit_pulse:
                unit_bursts = ax4.plot(np.where(np.array(p_solutions_unit[i]) > 0)[0], np.linspace(3*i+2,3*i+2,len(np.where(np.array(p_solutions_unit[i]) > 0)[0])), 'o', c='blue', alpha=.5, label = 'Predicted (Unit)')
                if i==0:
                    hand.append(unit_bursts[0])
            if cont:
                continuous_bursts = ax4.plot(np.where(np.array(p_solutions_cont[i]) > 0)[0], np.linspace(3*i+3,3*i+3,len(np.where(np.array(p_solutions_cont[i]) > 0)[0])), 'o', c='red', alpha=.5, label = 'Predicted (Continuous)')
                if i==0:
                    hand.append(continuous_bursts[0])
            if i == 0:
                a=2
                #plt.legend(handles=hand)

        plt.ylim([0,min(3*len(pulse_sizes_true)*1.1,150*1.1)])
        plt.grid(True)
        plt.xlabel('Time (Time Steps)')
        plt.ylabel('Trace Group')
        plt.title('Pulse Patterns Over Time')

        fig_pulse_scatter.savefig(os.path.join(outpath, subfolder_name, 'pulse_scatter.png'))
        #plt.show()

        #Pulse Size vs. freq (by trace)
        fig_size_freq_scatter = plt.figure(figsize=(15,5))
        ax4 = fig_size_freq_scatter.add_subplot(111)

        mean_size_true = []
        mean_size_unit = []
        mean_size_cont = []

        mean_freq_true = []
        mean_freq_unit = []
        mean_freq_cont = []

        for i in xrange(len(pulse_sizes_true)):
            mean_size_true.append(np.mean(np.array(pulse_sizes_true[i])))
            mean_freq_true.append(1.0 / (np.mean(np.array(on_lengths_true[i])) + np.mean(np.array(off_lengths_true[i]))))

            if unit_pulse:
                mean_size_unit.append(np.mean(np.array(pulse_sizes_unit[i])))
                mean_freq_unit.append(1.0 / (np.mean(np.array(on_lengths_unit[i])) + np.mean(np.array(off_lengths_unit[i]))))
            if cont:
                mean_size_unit.append(np.mean(np.array(pulse_sizes_cont[i])))
                mean_freq_unit.append(1.0 / (np.mean(np.array(on_lengths_cont[i])) + np.mean(np.array(off_lengths_cont[i]))))

        if simulation:
            true = ax4.plot(mean_freq_true, mean_size_true,'o',c='green', label = 'Actual')
        if unit_pulse:
            unit = ax4.plot(mean_freq_unit, mean_size_unit,'o',c='blue', label = 'Predicted (Unit)')
        if cont:
            continuous = ax4.plot(mean_freq_cont, mean_size_cont, 'o', c = 'red', label = 'Actual (Predicted (Continuous)')

        plt.grid(True)
        plt.xlabel('Mean Burst Frequency (1/(Time Step))')
        plt.ylabel('Mean Burst Size')
        plt.title('Mean Burst Size and Frequency by Trace')
        #plt.legend()
        plt.ylim([0,1.1*np.max(list(chain(*[mean_size_unit,mean_size_true])))])
        plt.xlim([0,1.1*np.max(list(chain(*[mean_freq_unit,mean_freq_true])))])
        fig_size_freq_scatter.savefig(os.path.join(outpath, subfolder_name, 'pulse_size_freq_scatter.png'))
        #plt.show()