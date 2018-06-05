% Script to generate Viterbi Fits for Inference traces
close all
clear 
addpath('E:\Nick\projects\hmmm\src\utilities\');
% addpath('D:\Data\Nick\projects\hmmm\src\utilities\');
%------------------------------Set System Params--------------------------%
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 30; % size of sliding window used
t_inf = 40;
%-----------------------------ID Variables--------------------------------%
w = 7; %memory assumed for inference
K = 3; %states used for final inference
Tres = 20; %time resolution
aggregate_fits = 0;  % if 1 apply same params to each trace regardless of stripe identity
%%% id variables
datatype = 'weka';
inference_type = 'dp';
project = 'eve7stripes_inf_2018_04_28'; %project identifier

%%% Generate filenames and writepath
id_string = [ '/w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/K' num2str(K) '_t_window' num2str(t_window) ...
     '_t_inf' num2str(round(t_inf)) '_' inference_type '/']; 

InfResultPath = ['../../dat/' project id_string];
OutPath = ['../../dat/' project id_string '/'];
mkdir(OutPath);
%%% Load Data
load([InfResultPath '\hmm_results_t_window' num2str(t_window) '_t_inf' num2str(t_inf) ...
    '.mat'])
load(['../../dat/' project '/inference_traces_' project '_dT' num2str(Tres) '.mat'])

inference_traces = trace_struct_final([trace_struct_final.inference_flag]==1);
particle_id_vec = [inference_traces.ParticleID];

%%% Viterbi Fits
hmm_regions = [hmm_results.binID];
alpha = hmm_results(1).alpha;
dT = hmm_results(1).dT;
viterbi_fit_struct = struct;
skipped_stripes = [];
i_pass = 1;

parfor i = 1:length(trace_struct_final)
%     MeanAP = round(trace_struct_final(i).MeanAP);        
    stripe_id =round(mode(trace_struct_final(i).stripe_id_vec),1);
    if aggregate_fits
        hmm_bin = hmm_results(round(hmm_regions,1)==0);
    else
%         hmm_index = find(MeanAP<=hmm_rb&MeanAP>=hmm_lb;        
        hmm_bin = hmm_results(hmm_regions==stripe_id);
    end        
    viterbi_fit_struct(i).skipped = 0;
    viterbi_fit_struct(i).hmm_bin = stripe_id;
    if isempty(hmm_bin)|| length(trace_struct_final(i).fluo_interp) < w        
%         error('wtf')
        viterbi_fit_struct(i).skipped = 1;
        viterbi_fit_struct(i).ParticleID = trace_struct_final(i).ParticleID;
        continue
    end    
    v = hmm_bin.initiation_mean/60*dT;
    noise = hmm_bin.noise_mean;
    pi0_log = log(ones(1,K)/3);
    A_log = reshape(log(hmm_bin.A_mean),K,K);                
    fluo = trace_struct_final(i).fluo_interp;            
    v_fit = viterbi (fluo, v', noise, pi0_log, A_log, K, w, alpha);        
    v_fit.time_exp = trace_struct_final(i).time_interp;
    v_fit.fluo_exp = fluo;            
    v_fit.v = v;
    v_fit.w = w;            
    v_fit.alpha = alpha;
    v_fit.noise = noise;
    v_fit.pi0 = exp(pi0_log);
    v_fit.A = exp(A_log);
    v_fit.agg_fit_flag = aggregate_fits;    
%     v_fit.trace_source = datapath;
    viterbi_fit_struct(i).v_fit = v_fit;
    viterbi_fit_struct(i).ParticleID = trace_struct_final(i).ParticleID;
%     i_pass = i_pass + 1;
    disp([num2str(i) ' of ' num2str(length(trace_struct_final)) ' completed'])
end
save([OutPath '/viterbi_fits_w' num2str(w) '_t_window' num2str(t_window) '_t_inf' num2str(t_inf) ...
    '.mat'],'viterbi_fit_struct') 