% Script to generate combined data set containing inference traces alongside
% inference results
addpath('../utilities/');
clear 
close all
%%%%%%-----Set System Params
w = 7; %memory assumed for inference
K = 2; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
fluo_field = 1;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 30;
t_inf = 35;
t_start_stripe = 10; 
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = [bin_range_vec stripe_range(i) + j/3 - 2/3];
    end
end

% id variables
datatype = 'weka';
inference_type = 'set';
project = 'eve7stripes_inf_2018_02_20'; %project identifier

%Generate filenames and writepath
% truncated_inference_w7_t20_alpha14_f1_cl1_no_ends1
id_thing = [ '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins)  '/']; 


DataPath = ['../../dat/' project '/' id_thing '/K' num2str(K) '_summary_stats/' ];
mkdir(DataPath)
% load inference summary 
load([DataPath 'hmm_results_mem' num2str(w) '_states' num2str(K) '_ti' num2str(t_inf)...
    '_tw' num2str(t_window) '.mat'])
% load viterbi fits
load([DataPath 'viterbi_fits_ti' num2str(t_inf) '_tw' num2str(t_window) '.mat'])
v_specific = viterbi_fit_struct;
load([DataPath 'viterbi_fits_ti25_tw50.mat'])
v_agg = viterbi_fit_struct;
% load traces
load('..\..\dat\eve7stripes_inf_2018_02_20\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat');

inference_traces = trace_struct_final([trace_struct_final.inference_flag]==1);
trace_particle_vec = [inference_traces.ParticleID];
inf_particle_index = [inference_traces.ParticleID];
inf_set_index = [inference_traces.setID];
va_particle_index = [v_agg.ParticleID];
vs_particle_index = [v_specific.ParticleID];
hmm_stripe_index = [hmm_results.binID];
for k = 1:length(inf_particle_index)            
    a = length(inf_particle_index) - k + 1;
    ParticleID = inf_particle_index(a);
    setID = inf_set_index(a);
    % extract relevant particle metrics    
    trace_time = inference_traces(inf_particle_index==ParticleID).time_interp;
    trace_fluo = inference_traces(inf_particle_index==ParticleID).fluo_interp;                
    stripe_id_vec = inference_traces(inf_particle_index==ParticleID).stripe_id_vec_interp(trace_time>t_start_stripe*60);
    if sum(isnan(stripe_id_vec)) > .5 * length(stripe_id_vec) || isempty(stripe_id_vec)
        continue
    end
    stripe_id = mode(stripe_id_vec);  
    hmm_bin = hmm_results(round(hmm_stripe_index,1)==round(stripe_id,1));
    R_mat = repmat(hmm_bin.R_fit_mean,length(trace_time),1);
    R_mat = R_mat(:,2:3)/60;
    r_mat = repmat(hmm_bin.initiation_mean'/60,length(trace_time),1);
    va_id = find(va_particle_index==ParticleID);
    if isempty(va_id)
        va_mat = NaN(length(trace_time),2);
    else
        v_fit_a = v_agg(va_id).v_fit;
        va_mat = [double(v_fit_a.z_viterbi') v_fit_a.fluo_viterbi'];
    end
    vs_id = find(vs_particle_index==ParticleID);
    if isempty(va_id)
        vs_mat = NaN(length(trace_time),2);
    else
        v_fit_s = v_specific(vs_id).v_fit;
        vs_mat = [double(v_fit_s.z_viterbi') v_fit_s.fluo_viterbi'];
    end
    
    pID_vec = repelem(ParticleID,length(trace_time))';
    setID_vec = repelem(setID,length(trace_time))';
    stripe_id_long = repelem(stripe_id,length(trace_time))';
    
    trace_mat = [round(pID_vec,4) double(setID_vec) double(stripe_id_long) double(trace_time') trace_fluo' ...
        double(vs_mat) double(va_mat) R_mat double(r_mat)];
    trace_mat = trace_mat(1:end-w-1,:);        
    if k == 1
        trace_mat_long = trace_mat;
    else
        trace_mat_long = [trace_mat_long ; trace_mat];
    end
end
% remove NaN rows
% trace_mat_long = trace_mat_long(nan_index,:);
header = {'particle_id', 'set_id', 'stripe_id','time', 'fluo', 'v_state_sp', 'v_fluo_sp'...
          'v_state_agg', 'v_fluo_agg', 'k_on', 'k_off', 'initiation_rate_off', 'initiation_rate_on'};
csvwrite_with_headers([DataPath '\eve_data_longform.csv'], ...
                       trace_mat_long, header,9); 
                 
save([DataPath 'eve_data_longform.mat'],'trace_mat_long')                   