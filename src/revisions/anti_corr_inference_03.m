% Script to Conduct Windowed HMM Inference Validation Tests
% Now explicitly focused high rate-induced artifacts
close all
clear 
addpath('../../../hmmm/src/utilities');
% Specify whether inference is being conducted on Savio Cluster
savio = 0; 
if savio
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    rep_num_list = zeros(1,length(savio_groups));
    for i =1:length(rep_num_list)
        rep_num_list(i) = savio_groups{i};
    end
    disp(rep_num_list);
else
    rep_num_list = 1:10;
end
warning('off','all') %Shut off Warnings
% Load simulation data
DataPath = '../../dat/revisions/anti_corr/';
load([DataPath 'sim_data_03.mat'])
FigPath = '../../fig/revisions/anti_corr/';
mkdir(FigPath);
%------------------Define Inference Variables------------------------------%
if savio
    MaxWorkers = 24;
else
    MaxWorkers = 12;
end
n_localEM = 25; % set num local runs
n_steps_max = 500; % set max steps per inference
eps = 10e-4; % set convergence criteria

K = sim_struct(1).K; % State(s) to use for inference
w = sim_struct(1).w; % system memory

%-------------------------------System Vars-------------------------------%
dT = sim_struct(1).dT; % Time Resolution
alpha = sim_struct(1).alpha;
inference_times = 20*60;
t_window = 20*60; % determines width of sliding window
project = 'anti_corr_03/';
%---------------------- Set write path -----------------------------------%
out_suffix =  ['/' project '/w' num2str(w) '/states' num2str(K) '/']; 
if savio
    out_prefix = '/global/scratch/nlammers/';
else    
    out_prefix = '../../out/revisions/';
end
%%%Set write directories
out_dir = [out_prefix out_suffix];
mkdir(out_dir);

%% Conduct Inference
% structure array to store the analysis data
n_groups = length(sim_struct);
local_meta = struct; % Store local_EM results
init_meta = struct; % Store initiation info
n_bootstraps = 5;
boot_size = 100; % n_traces
for g = 1:n_groups % loop through different Inference groups
    sim_name = sim_struct(g).ID;
    fluo_data_full = sim_struct(g).fluo_data;                
    for t = 1:length(inference_times)
        t_inf = inference_times(t);
        t_start = max(t_inf - t_window/2,(w+1)*dT);
        t_stop = t_inf + t_window/2;
        trace_time = sim_struct(g).trace_time;
        index_vec = 1:length(trace_time);
        t_filter = trace_time>=t_start&trace_time<t_stop;            
        for b = 1:n_bootstraps% multiple versions of each replicate                
            iter_start = now;            
            local_struct = struct;
            init_struct = struct;
            output = struct;            
            % Use current time as unique inference identifier 
            inference_id = num2str(round(10e5*now));            
            % Generate filenames            
            fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) ...
                '_' sim_name 'boot_' num2str(b) '_' sim_struct(g).ID ...
                '_time' num2str(round(t_inf/60)) '_id' inference_id];                
            out_file = [out_dir '/' fName_sub];
            % Initialize logL to - infinity
            logL_max = -Inf;
            % Extract fluo_data                
            fluo_data = cell([boot_size, 1]);                
            boot_indices = randsample(1:length(fluo_data_full),boot_size,true);
            iter = 1;
            for f = boot_indices
                f_trace = fluo_data_full{f};                                    
                fluo_data{iter} = f_trace;                    
                iter = iter + 1;
            end                                
            n_dp = length([fluo_data{:}]);                 
            % Random initialization of model parameters
            param_init = initialize_random (K, w, fluo_data);
            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory_truncated (fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, 500, 1e-4);                
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);            
            p = gcp('nocreate');
            if isempty(p)
                parpool(MaxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > MaxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(MaxWorkers);
            end
            parfor i_local = 1:n_localEM % Parallel Local EM                
                % Random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;
                noise_init = param_init.noise;
                % Record
                init_struct(i_local).A_init = exp(A_log_init);                
                init_struct(i_local).v_init = v_init;
                init_struct(i_local).noise_init = noise_init;                
                init_struct(i_local).subset_id = i_local;
                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);                    
                %---------------------------------------------------------%                
                % Save Results 
                local_struct(i_local).inference_id = inference_id;
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;                
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / dT;                                
                local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct(i_local).pi0 = exp(local_out.pi0_log);
%                         local_struct(i_local).total_time = local_out.runtime;
                local_struct(i_local).total_steps = local_out.n_iter;               
                local_struct(i_local).soft_struct = local_out.soft_struct;               
            end
            [logL, max_index] = max([local_struct.logL]); % Get index of best result           
%                 output.local_runs = local_struct;            
            % Save parameters from most likely local run
            output.soft_struct = local_struct(max_index).soft_struct;                        
            output.pi0 = local_struct(max_index).pi0;                                        
            output.r = local_struct(max_index).r(:);            
            output.noise = local_struct(max_index).noise;
            output.A = local_struct(max_index).A(:);
            output.A_mat = local_struct(max_index).A;                                    
            % Info about run time
            output.total_steps = local_struct(max_index).total_steps;                                  
            output.total_time = 100000*(now - iter_start);            
            % Save inference ID variables
            output.sim_name = sim_name;                            
            output.boot_num = b;
            output.boot_indices = boot_indices;
            % time info
            output.t_window = t_window;
            output.t_inf = t_inf;                            
            % Other inference characteristics                            
            output.n_dp = n_dp;
            output.w = w;
            output.alpha = alpha;
            output.deltaT = dT;                    
            save([out_file '.mat'], 'output');           
        end
    end  
end