% Script to generate Viterbi Fits for Inference traces
% addpath('D:\Data\Nick\projects\hmmm\src\utilities\');
addpath('E:\Nick\projects\hmmm\src\utilities\');
%------------------------------Set System Params--------------------------%
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 50; % size of sliding window used
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

%%% Generate filenames and writepath
id_string = [ '/w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/states' num2str(K) '_t_window' num2str(t_window)  '_' inference_type '/']; 

InfResultPath = ['../../dat/' project '/' id_string];
%%% Load Data
load('..\..\dat\eve7stripes_inf_2018_02_20\w7_t20_alpha14_f1_cl1_no_ends1_tbins1\states2\t_window30\set_bootstrap_results\hmm_results_mem7_states2.mat') %Inference results
load('..\..\dat\eve7stripes_inf_2018_02_20\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat')
%%
% load(['../../dat/' project '/inference_traces_' project '.mat']) % load traces
inference_traces = trace_struct_final([trace_struct_final.inference_flag]==1);
particle_id_vec = [inference_traces.ParticleID];

%%% Viterbi Fits
hmm_regions = [hmm_results.binID];
alpha = hmm_results(1).alpha;
dT = hmm_results(1).dT;
viterbi_fit_struct = struct;
skipped_stripes = [];
for i = 1:length(inference_traces)
    stripe_id_vec = inference_traces(i).stripe_id_vec_interp;
    stripe_id = round(mode(stripe_id_vec),1);   
    hmm_bin = hmm_results(round(hmm_regions,1)==stripe_id);
    if isempty(hmm_bin)
        skipped_stripes = [skipped_stripes stripe_id];
        v_fit = struct;
        v_fit.ParticleID = inference_traces(i).ParticleID; 
        viterbi_fit_struct(i).ParticleID = inference_traces(i).ParticleID; 
%         warning('No inference results found for region. Skipping')
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
%     v_fit.trace_source = datapath;
    viterbi_fit_struct(i).v_fit = v_fit;
    viterbi_fit_struct(i).ParticleID = inference_traces(i).ParticleID;
    disp([num2str(i) ' of ' num2str(length(inference_traces)) ' completed'])
end

save([InfResultPath '/viterbi_fits.mat'],'viterbi_fit_struct') 