from matplotlib import pyplot as plt
from scipy.optimize import least_squares
import numpy as np
from pulp import *
from itertools import chain
from generate_traces import generate_traces_gill, generate_traces_prob_mat

def shape_mode_decomposition(memory, length, input_size, batch_size, out_mem=1, v_num=None,r_mat=np.array([]), v=np.array([]), switch_low=1, switch_high = 12, noise_scale =.025, alpha=1.0):

    #Generate set of traces using specified rate matrix and emission vector
    _, _, _, _, fluo_list, _ =  generate_traces_gill(memory, length, input_size, batch_size,
                                                     num_steps=1,
                                                     r_mat=r_mat,
                                                     v=v,
                                                     noise_scale =.025,
                                                     alpha=1.0)

    #Break traces into fragments of length w and store in list
    fragments = []

    for i in xrange(batch_size):
        n_iters = length / memory
        for j in xrange(n_iters-1):
            fragments.append(np.round(fluo_list[j],3)[j*memory:(j+1)*memory].tolist())

    #Convert list to array
    fragment_array = np.array(fragments)
    #Subtract first element

    fragment_array = fragment_array - np.matrix.transpose(np.tile(fragment_array[:,0],(memory,1)))

    #Subtract mean
    fragment_mean = np.mean(fragment_array,axis=0)
    frag_diff_array = fragment_array - np.tile(fragment_mean,(len(fragment_array[:,0]),1))
    #Get covariance matrix
    fragment_correlations = np.matmul(np.matrix.transpose(frag_diff_array),frag_diff_array)

    #Get eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eig(fragment_correlations)

    return(eigenvalues, eigenvectors, fragment_mean)

def lsq_fit(fragment, shape_mat, mean_shape):
    def loss_fun(x,fragment,shape_mat, mean_shape):
        return np.sum(np.square(fragment - np.sum(np.matrix.transpose(np.tile(x,(len(shape_mat[:,0]),1)))*np.matrix.transpose(shape_mat),axis=0) - mean_shape))

    # Define initialization vector for states
    x0 = np.ones(len(shape_mat[:,0]))
    res = least_squares(loss_fun, x0, args=(fragment, shape_mat, mean_shape), verbose=1)
    return res.x

def shape_mode_fits(shape_modes, mean_shape, traces, n_modes):
    #Takes 2D array of raw traces and returns 3D array with third D corresponding to eigenvectors
    #Traces are fit piecewise in segments equal to length of system memory
    memory = len(shape_modes)
    n_traces = len(traces[0,:])
    mem_diff = len(traces[:,0]) - (len(traces[:,0]) / memory) * memory
    #Let's pad trace tails with zeros to ensure they
    trace_mat = np.concatenate((traces,np.zeros((mem_diff,len(traces[0,:])))))
    t = len(trace_mat[:,0])
    #Generate 3D output array with third dimension equal to num eigenvectors + mean + basis value
    trace_mat_out = np.zeros((t,n_traces,memory+2))
    tr_fits = []
    #Now we fit traces piecewise
    for i in xrange(n_traces):
        n_iters = t / memory
        trace_fit = np.zeros((t,memory+2))
        for j in xrange(n_iters):
            weight_vec = lsq_fit(fragment=trace_mat[j*memory:(j+1)*memory,i] - trace_mat[j*memory,i],
                                 shape_mat=shape_modes,
                                 mean_shape=mean_shape)
            #print(weight_vec)
            trace_fit[j*memory:(j+1)*memory,0] = trace_mat[j*memory,i]
            trace_fit[j*memory:(j+1)*memory,1] = mean_shape
            trace_fit[j*memory:(j+1)*memory,2:] = shape_modes*np.tile(weight_vec,(memory,1))
        trace_mat_out[:,i,:] = trace_fit
        tr_fits.append(np.sum(trace_fit[0:len(traces[:,0]),0:n_modes],axis=1).tolist())
    return trace_mat_out, tr_fits

#Look for consecutive 0 and nonzero elements in promoter state vector
def on_off_streaks(data):
    z_lengths = []
    nz_lengths = []
    ind_out = []
    if len(data) > 0:
        seq_len = len(data)
        bin_vec = np.array(1*(data > 0))
        bin_diff_vec = np.diff(bin_vec)
        ind_vec = np.where(bin_diff_vec != 0)[0] + 1
        streak_vec = ind_vec - np.array(list(chain(*[[0],ind_vec[:-1]])))
        streak_vec = streak_vec.tolist()
        if seq_len != ind_vec[-1]:
            streak_vec.append(seq_len - ind_vec[-1])
        nz_id = 0
        z_id = 0
        if bin_diff_vec[ind_vec[0]-1] == 1:
            nz_id = 1
        elif bin_diff_vec[ind_vec[0]-1] == -1:
            z_id = 1

        z_lengths = [streak_vec[z_id + 2*i] for i in xrange(len(streak_vec)/2 + (len(streak_vec) % 2)*(1 - z_id))]
        nz_lengths = [streak_vec[nz_id + 2*i] for i in xrange(len(streak_vec)/2 + (len(streak_vec) % 2)*(1 - nz_id))]
        ind_out =  np.array(list(chain(*[[0],ind_vec,[seq_len+1]])))

    return(z_lengths, nz_lengths,ind_out)

def get_pulse_size(data, change_points):
    pulse_size_list = []
    for i in xrange(len(change_points)-1):
        pulse_size_list.append(np.sum(data[change_points[i]:change_points[i+1]]))
    pulse_size_list = np.array(pulse_size_list)
    pulse_size_list = pulse_size_list[pulse_size_list!=0].tolist()

    return(pulse_size_list)


if __name__ == "__main__":
    """
    # memory
    w = 15
    # Fix trace length for now
    T = 360
    # Input magnitude
    F = 501
    # Number of traces per batch
    batch_size = 250
    R = np.array([[-.007,.007,.006],[.004,-.01,.008],[.003,.003,-.014]]) * 6.0
    v = np.array([0.0,4.0,8.0])
    eigenvalues, eigenvectors, fragment_mean = shape_mode_decomposition(w,length=T, input_size=w*max(v)+1,batch_size=batch_size, r_mat=R, v=v, noise_scale =.025, alpha=1.0)
    #print(np.real(eigenvalues))
    #sys.exit(1)
    _, _, _, _, fluo_list, _ = generate_traces_gill(w,length=T, input_size=w*max(v)+1,batch_size=1, r_mat=R, v=v, noise_scale =.025, alpha=1.0,num_steps=1)
    #print(eigenvectors[0])
    traces = np.array(fluo_list)
    test, fit = shape_mode_fits(np.real(eigenvectors), fragment_mean, np.matrix.transpose(traces), n_modes)
    slice = np.sum(test[:,0,0:6],axis=1)

    fig_fluo = plt.figure(figsize=(12, 4))
    ax = plt.subplot(1, 1, 1)
    ax.plot(slice, 'b-o', alpha=0.4, label='Actual')
    #print(fit[0])
    ax.plot(traces[0], c='r', label='Predicted')
    #ax.plot(np.array(fit[0]),'ro' , label='Actual')
    plt.show()
    """
    vec = np.array([0,0,0,0,23,354,23,4,0,0,2,3,0])

    z, nz, change_points = on_off_streaks(vec)
    pulse_sizes = get_pulse_size(vec, change_points)
    print(z)
    print(nz)
    print(change_points)
    print(pulse_sizes)
