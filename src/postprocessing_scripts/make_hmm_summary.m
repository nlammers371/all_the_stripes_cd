% Script to make HMM parameter csv files
addpath('../utilities/');
clear 
close all
%%%%%%-----Set System Params
w = 7; %memory assumed for inference
K = 3; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
fluo_field = 1;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 30;
t_inf = 40;
t_start_stripe = 25; 
%-----------------------------ID Variables--------------------------------%

% id variables
datatype = 'weka';
inference_type = 'dp';
project = 'eve7stripes_inf_2018_04_28'; %project identifier

%Generate filenames and writepath
id_thing = [ '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins)  '/']; 


DataPath = ['../../dat/' project '/' id_thing '/K' num2str(K) '_summary_stats/' ];
% load inference summary 
load([DataPath 'hmm_results_t_window' num2str(t_window) '_t_inf' num2str(t_inf)  '.mat'])
%%
hmm_stripe_index = [hmm_results.binID]; % inf region ref vec

header = {'stripe_id','k_on', 'k_off', 'initiation_rate_off', 'initiation_rate_on'};
param_mat = NaN(length(hmm_stripe_index), length(header));
for i = 1:length(hmm_stripe_index)
    param_mat(i,1) = hmm_stripe_index(i);
    R_full = hmm_results(i).R_fit_mean;
    if K == 2
        R = R_full(2:3);
        r = hmm_results(i).initiation_mean';
    elseif K == 3
        effective_off = ((R_full(4)+R_full(6))./R_full(4).*-R_full(5).^-1 + ...
                        R_full(6) ./ R_full(4).*-R_full(9).^-1).^-1;
        R = [R_full(2) effective_off];                
        r_full = hmm_results(i).initiation_mean;
        s2 = hmm_results(i).occupancy_mean(3);
        s1 = hmm_results(i).occupancy_mean(2);
        eff_init = r_full(2);%+ s2/(s1+s2)*r_full(3);
        r = [r_full(1) eff_init];                
    end
    param_mat(i,2:end) = [R r];
end
%%
csvwrite_with_headers([DataPath '\eve_hmm_summary_K' num2str(K) '.csv'], ...
                       param_mat, header,9); 
                 
save([DataPath 'eve_hmm_summary_K' num2str(K) '.mat'],'param_mat') 