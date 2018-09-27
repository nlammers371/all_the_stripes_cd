% script to generate data sets for testing possible link between detection
% threshold and apparent transition rate correlation
clear 
close all
% set paths
addpath('../../../hmmm/src/utilities');
ResultsPath = '..\..\dat\_initial_eLife_submission\eve7stripes_inf_2018_04_28\w7_t20_alpha14_f1_cl1_no_ends1_tbins1\K3_summary_stats\';
DataPath = '../../dat/revisions/anti_corr/';
mkdir(DataPath);
load([ResultsPath 'hmm_results_t_window30_t_inf40.mat'])
FigPath = '../../fig/revisions/anti_corr/';
mkdir(FigPath);
% extract average bursting parameters
bin_index = [hmm_results.binID];
K = 3; % N states
w = 5; % memory
R_mean = reshape(hmm_results(bin_index==0).R_fit_mean,K,K)/60; % transition rates
r_mean = hmm_results(bin_index==0).initiation_mean/60; % initiation rates
r_mean = [0 r_mean(2) 2*r_mean(2)];
pi0 = [0 1 1]/2; % initial state PDF
alpha = w/7*hmm_results(bin_index==0).alpha; %MS2 rise time
dT = hmm_results(bin_index==0).dT; % time resolution
noise_mean = w/7*hmm_results(1).noise_mean; % noise level (AU)
seq_length = 60/dT*40;
trace_time = (1:seq_length)*20;
% set simulation parameters
n_traces = 1000;
granularity = 100;
% Use average inference results to graded expression using only kon
% or only koff. For simplicity assume independent promoters for data
% generation
k_on_bounds = [.005 .015];
k_off_bounds = [.005 .015];
sim_ids = {'random','anticorrelated'};
sim_struct = struct;
% generate sample traces
for i = 1:2    
    traces = cell(n_traces,1);
    trajectories = cell(n_traces,1);
    jumps = cell(n_traces,1);
    f_vec = NaN(1,n_traces);
    k_on_vec = NaN(1,n_traces);
    k_off_vec = NaN(1,n_traces);
    for n = 1:n_traces     
        if i == 1
            kon = rand()*(k_on_bounds(2)-k_on_bounds(1))+k_on_bounds(1);
            koff = rand()*(k_off_bounds(2)-k_off_bounds(1))+k_off_bounds(1);
        elseif i == 2
            kon = rand()*(k_on_bounds(2)-k_on_bounds(1))+k_on_bounds(1);
            koff = (1-(kon-k_on_bounds(1))/(k_on_bounds(2)-k_on_bounds(1))) ...
                * (k_off_bounds(2)-k_off_bounds(1)) + k_off_bounds(1);
        end        
        R = [-2*kon koff 0 ; 2*kon -koff-kon 2*koff; 0 kon -2*koff];            
        gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                        K, w, R, dT, r_mean, noise_mean, pi0);
        traces{n} = gillespie.fluo_MS2;
        trajectories{n} = gillespie.naive_states;
        jumps{n} = gillespie.transition_times;
        k_on_vec(n) = kon;
        k_off_vec(n) = koff;
        f_vec(n) = nanmean(traces{n});        
    end
    sim_struct(i).ID = sim_ids{i};    
    sim_struct(i).dT = dT;
    sim_struct(i).fluo_data  = traces;
    sim_struct(i).trace_time = trace_time;
    sim_struct(i).promoter_states  = trajectories;
    sim_struct(i).jump_times = jumps;    
    sim_struct(i).r_emission = r_mean;
    sim_struct(i).noise = noise_mean;
    sim_struct(i).w = w;
    sim_struct(i).alpha = alpha;
    sim_struct(i).K = K;
    sim_struct(i).k_on_bounds = k_on_bounds;
    sim_struct(i).k_off_bounds = k_off_bounds;
    sim_struct(i).k_on_bounds = k_on_bounds;
    sim_struct(i).k_off_vec = k_off_vec;
    sim_struct(i).k_on_vec = k_on_vec;
    
    rate_scatter = figure;    
    scatter(k_on_vec,k_off_vec,40,f_vec,'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
    xlabel('k_{on} (s^{-1})')
    ylabel('k_{off} (s^{-1})')
    title(['k_{on} vs. k_{off} (' sim_ids{i} ')'])
    colorbar
    saveas(rate_scatter,[FigPath 'rate_scatter_' sim_ids{i} '.png'])
        
end
save([DataPath 'sim_data_01.mat'],'sim_struct')