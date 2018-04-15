% Script to generate Viterbi Fits for Inference traces

close all
clear 
addpath('E:\Nick\projects\hmmm\src\utilities\');
% addpath('D:\Data\Nick\projects\hmmm\src\utilities\');
%------------------------------Set System Params--------------------------%
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 30; % size of sliding window used
t_inf = 35;
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = stripe_range(i) - j/3 + 2/3;
    end
end


w = 7; %memory assumed for inference
K = 2; %states used for final inference
Tres = 20; %time resolution

% id variables
datatype = 'weka';
inference_type = 'dp';
project = 'eve7stripes_inf_2018_02_20'; %project identifier
aggregate_fits = 1; 
%%% Generate filenames and writepath
id_string = [ '/w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/K' num2str(K) '_t_window' num2str(t_window)  '_' inference_type '/']; 

InfResultPath = ['../../dat/' project id_string];
OutPath = ['../../dat/' project id_string '/K' num2str(K) '_viterbi_fits/'];
mkdir(OutPath);
%%% Load Data
% load([InfResultPath '\hmm_results_mem' num2str(w) '_states' num2str(2) '.mat']) %Inference results
load('E:\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_02_20\w7_t20_alpha14_f1_cl1_no_ends1_tbins1\K2_summary_stats\hmm_results_mem7_states2_ti25_tw50.mat')
load('..\..\dat\eve7stripes_inf_2018_02_20\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat')

% load(['../../dat/' project '/inference_traces_' project '.mat']) % load traces
inference_traces = trace_struct_final([trace_struct_final.inference_flag]==1);
particle_id_vec = [inference_traces.ParticleID];

%%% Viterbi Fits
hmm_regions = [hmm_results.binID];
alpha = hmm_results(1).alpha;
dT = hmm_results(1).dT;
viterbi_fit_struct = struct;
skipped_stripes = [];
i_pass = 1;
for i = 1:length(inference_traces)
    stripe_id_vec = inference_traces(i).stripe_id_vec_interp;
    stripe_id = round(mode(stripe_id_vec),1);
    if aggregate_fits
        hmm_bin = hmm_results(round(hmm_regions,1)==0);
    else
        hmm_bin = hmm_results(round(hmm_regions,1)==stripe_id);
    end
    if isempty(hmm_bin)        
        continue
    end    
    v = hmm_bin.initiation_mean/60*dT;
    noise = hmm_bin.noise_mean;
    pi0_log = log([.5;.5]);
    A_log = reshape(log(hmm_bin.A_mean),K,K);                
    fluo = inference_traces(i).fluo_interp;            
    v_fit = viterbi (fluo, v', noise, pi0_log, ...
                                    A_log, K, w, alpha);    
    v_fit.time_exp = inference_traces(i).time_interp;
    v_fit.fluo_exp = fluo;            
    v_fit.v = v;
    v_fit.w = w;        
    v_fit.alpha = alpha;
    v_fit.noise = noise;
    v_fit.pi0 = exp(pi0_log);
    v_fit.A = exp(A_log);
    v_fit.agg_fit_flag = aggregate_fits;
%     error('afsa')
%     v_fit.trace_source = datapath;
    viterbi_fit_struct(i_pass).v_fit = v_fit;
    viterbi_fit_struct(i_pass).ParticleID = inference_traces(i).ParticleID;
    i_pass = i_pass + 1;
    disp([num2str(i) ' of ' num2str(length(inference_traces)) ' completed'])
end
save([OutPath '/viterbi_fits_ti' num2str(t_inf) '_tw' num2str(t_window)...
    '.mat'],'viterbi_fit_struct') 