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
% set simulation parameters
n_traces = 1000;
granularity = 100;
% Use average inference results to graded expression using only kon
% or only koff. For simplicity assume independent promoters for data
% generation
k_on_mean = .5*(R_mean(2,1)/2 + R_mean(3,2));
k_off_mean = .5*(R_mean(1,2) + R_mean(2,3)/2);
%%
%%% kon %%%
k_on_min = 3/4*k_on_mean;
k_on_max = 3/2*k_on_mean;
k_on_vec = linspace(k_on_min,k_on_max,100);
ap_vec = 1:100;
% sample from k_on_vec
sim_struct = struct;
k_on_ap = NaN(1,n_traces);
k_on_traces = cell(1,n_traces);
k_on_trajectories = cell(1,n_traces);
k_on_jumps = cell(1,n_traces);
k_on_f_vec = NaN(1,n_traces);
for n = 1:n_traces
    ap = randsample(ap_vec,1);
    kon = k_on_vec(ap);
    R = [-2*kon k_off_mean 0 ; 2*kon -k_off_mean-kon 2*k_off_mean; 0 kon -2*k_off_mean];            
    gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                    K, w, R, dT, r_mean, noise_mean, pi0);
    k_on_traces{n} = gillespie.fluo_MS2;
    k_on_trajectories{n} = gillespie.naive_states;
    k_on_jumps{n} = gillespie.transition_times;
    k_on_ap(n) = ap;
    k_on_f_vec(n) = nanmean(k_on_traces{n});
end
sim_struct(1).ID = 'k_on_traces';
sim_struct(1).dT = dT;
sim_struct(1).fluo_data  = k_on_traces;
sim_struct(1).promoter_states  = k_on_trajectories;
sim_struct(1).jump_times = k_on_jumps;
sim_struct(1).ap_vec = k_on_ap;
sim_struct(1).r_emission = r_mean;
sim_struct(1).noise = noise_mean;
sim_struct(1).w = w;
sim_struct(1).alpha = alpha;
sim_struct(1).K = K;
sim_struct(1).k_on_mean = k_on_mean;
sim_struct(1).k_off_mean = k_off_mean;
sim_struct(1).k_vec = k_on_vec;

%%% kon %%%
k_off_min = (2*k_off_mean-k_on_mean)/3;
k_off_max = k_on_mean + 2*k_off_min;
k_off_vec = linspace(k_off_min,k_off_max,100);
% sample from k_on_vec
k_off_ap = NaN(1,n_traces);
k_off_traces = cell(1,n_traces);
k_off_trajectories = cell(1,n_traces);
k_off_jumps = cell(1,n_traces);
k_off_f_vec = NaN(1,n_traces);
for n = 1:n_traces
    ap = randsample(ap_vec,1);
    koff = k_off_vec(ap);
    R = [-2*k_on_mean koff 0 ; 2*k_on_mean -koff-k_on_mean 2*koff; 0 k_on_mean -2*koff];            
    gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                    K, w, R, dT, r_mean, noise_mean, pi0);
    k_off_traces{n} = gillespie.fluo_MS2;
    k_off_trajectories{n} = gillespie.naive_states;
    k_off_jumps{n} = gillespie.transition_times;
    k_off_ap(n) = ap;
    k_off_f_vec(n) = nanmean(k_off_traces{n});
end
sim_struct(2).ID = 'k_off_traces';
sim_struct(2).dT = dT;
sim_struct(2).fluo_data  = k_off_traces;
sim_struct(2).promoter_states  = k_off_trajectories;
sim_struct(2).jump_times = k_off_jumps;
sim_struct(2).ap_vec = 100-k_off_ap+1;
sim_struct(2).r_emission = r_mean;
sim_struct(2).noise = noise_mean;
sim_struct(2).w = w;
sim_struct(2).alpha = alpha;
sim_struct(2).K = K;
sim_struct(2).k_on_mean = k_on_mean;
sim_struct(2).k_off_mean = k_off_mean;
sim_struct(2).k_vec = k_off_vec;

%% Make summary figure
cm = jet(128);
kon_scatter = figure;
hold on
scatter(k_on_vec(k_on_ap),k_on_f_vec,40,'MarkerFaceColor',cm(30,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
title('Trace Fluorescence vs. k_{on}')
xlabel('k_{on} (s^{-1})')
ylabel('AU')
saveas(kon_scatter,[FigPath 'k_on_fluo.png'])

koff_scatter = figure;
hold on
scatter(k_off_vec(k_off_ap),k_off_f_vec,40,'MarkerFaceColor',cm(110,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
title('Trace Fluorescence vs. k_{off}')
xlabel('k_{off} (s^{-1})')
ylabel('AU')
saveas(koff_scatter,[FigPath 'k_off_fluo.png'])

save([DataPath 'sim_data_01.mat'],'sim_struct');
close all