% Script to extract average parameters from HMM inference and use to make
% Viterbi Fits
close all
clear 
addpath('../../../hmmm/src/utilities');

% Load simulation data
DataPath = '../../dat/revisions/anti_corr/';
load([DataPath 'sim_data_01.mat'])
% Figure write path
FigPath = '../../fig/revisions/anti_corr/';
mkdir(FigPath);
% Load inference results
K = sim_struct(1).K; % N States
w = sim_struct(1).w; % Memory
dT = sim_struct(1).dT; % Time Resolution
alpha = sim_struct(1).alpha;
inference_times = 20*60;
t_window = 20*60; % determines width of sliding window
project = 'anti_corr_01/';
% Set read path
out_suffix =  ['/' project '/w' num2str(w) '/states' num2str(K) '/']; 
out_prefix = '../../out/revisions/';
out_dir = [out_prefix out_suffix];
%%
% Get list of inference files
file_list = dir([out_dir '*.mat']);
inference_results = struct;
rand_vec = [];
for f = 1:numel(file_list)
    load([out_dir file_list(f).name])
    fn = fieldnames(output);
    for i = 1:numel(fn)
        inference_results(f).(fn{i}) = output.(fn{i});
    end
    rand_vec = [rand_vec strcmp(output.sim_name,'random')];
end
rand_index = unique(rand_vec);
sigma_vec = [inference_results.detection_threshold];
sigma_index = unique(sigma_vec);
hmm_results = struct;
for i = 1:numel(sigma_index)    
    for k = 1:numel(rand_index)
        ind = (i-1)*numel(rand_index) + k;
        s_results = inference_results(rand_vec==rand_index(k)&sigma_vec==sigma_index(i));
        A_array = NaN(K,K,numel(s_results));
        v_array = NaN(K,numel(s_results));
        for j = 1:numel(s_results)
            v_vec = s_results(j).r*dT;
            [v_vec, vi] = sort(v_vec);
            A_mat = s_results(j).A_mat;
            A_array(:,:,j) = A_mat(vi,vi);
            v_array(:,j) = v_vec;
        end
        hmm_results(ind).v_mean = nanmean(v_array,2);
        hmm_results(ind).v_std = nanstd(v_array,[],2);
        A_mean = nanmean(A_array,3);
        A_std = nanstd(A_array,[],3);
        A_mean = A_mean ./ sum(A_mean);
        A_std = A_std ./ sum(A_mean);
        sigma_mean = nanmean(sqrt([s_results.noise]));
        % store
        hmm_results(ind).A_mean = A_mean;
        hmm_results(ind).A_std = A_std;        
        hmm_results(ind).w = w;
        hmm_results(ind).K = K;
        hmm_results(ind).dT = dT;
        hmm_results(ind).random = rand_index(k);
        hmm_results(ind).sigma_thresh = sigma_index(i);
        hmm_results(ind).noise_mean = sigma_mean;        
    end
end

% Save
save([DataPath 'sim_inf_results_w' num2str(w) '_K' num2str(K) '.mat'],'hmm_results')
hmm_sigma = [hmm_results.sigma_thresh];
hmm_random = [hmm_results.random];

% Now perform Viterbi Fitting
viterbi_fit_struct = struct;

for i = 1:length(sigma_index) % iterate through different thresholds
    sigma_thresh = sigma_index(i);
    for j = 1:length(rand_index) % iterate through diff sim types
        random = rand_index(j);
        fluo_data = sim_struct(j).fluo_data;
        hmm_bin = hmm_results(hmm_random==rand_index(j)&hmm_sigma==sigma_thresh);        
        ind_base = (i-1)*numel(rand_index)+j;
        v_fit = struct;
        parfor k = 1:numel(fluo_data) % iterate through traces            
            fluo = fluo_data{k};
            fluo(fluo<=sigma_thresh) = 0;    
            v_fit(k).random = rand_index(j);
            v_fit(k).sigma_thresh = sigma_thresh;    
            v_fit(k).trace_num = k;            
            v = hmm_bin.v_mean;
            noise = hmm_bin.noise_mean;
            pi0_log = log(ones(1,K)/3);
            A_log = reshape(log(hmm_bin.A_mean),K,K);                            
            vf = viterbi (fluo, v', noise, pi0_log, A_log, K, w, alpha);                    
            vf.time_exp = sim_struct(j).trace_time;
            vf.fluo_thresh = fluo;            
            vf.fluo_orig = fluo_data{k};            
            vf.v = v;
            vf.w = w;            
            vf.alpha = alpha;
            vf.noise = noise;
            vf.pi0 = exp(pi0_log);
            vf.A = exp(A_log);            
            fn = fieldnames(vf);
            for n = 1:numel(fn)
                v_fit(k).(fn{n}) = vf.(fn{n});
            end
            v_fit(k).sigma_thresh = sigma_thresh;
            v_fit(k).random = random;                        
        end
        viterbi_fit_struct(ind_base).sigma_thresh = sigma_thresh;
        viterbi_fit_struct(ind_base).random = random;
        viterbi_fit_struct(ind_base).v_fits = v_fit;
        disp(['Completed Iteration ' num2str(ind_base) ' of ' num2str(numel(rand_index)*numel(sigma_index))])
    end   
end

save([DataPath '/viterbi_fits_w' num2str(w) '_K' num2str(K) ...
    '.mat'],'viterbi_fit_struct') 