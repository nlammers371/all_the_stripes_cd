%Get environment variable from job script
% stripe_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};
stripe_groups = {1,2,3,4,5,6,7};
stripe_regions = [0];
%only do stripe centers for now
w = 7;
state_vec = [2];
%route to utilities folder
addpath('../utilities');
%------------------Define Inference Variables------------------------------%
% max num workers
pool_max = 12;
% set num local runs
n_localEM = 25;
%Time Resolution
Tres = 20;
% set max steps per inference
n_steps_max = 500;
% set convergence criteria
eps = 10e-4;
%Project ID
project = 'eve7stripes_inf_2017_09_25';

%Path to raw data
datapath = ['../../dat/' project '/'];
start_time = 0;
stop_time = 60;
% generate save names
dataname = ['inference_traces_t' num2str(Tres) '_' project '.mat'];

date_str = '2017-09-25';
out_dir_single =  ['../../out/' project '/' date_str '/' 'inference_w' ...
    num2str(w)  'set_boots/'];

if exist(out_dir_single) ~= 7
    mkdir(out_dir_single);
end

% Load data for inference into struct named: interp_struct
load([datapath dataname]);

%Get Tres and alpha
alpha = interp_struct(1).alpha_frac*w;
deltaT = Tres;

% create position index
stripeIndex = unique([interp_struct.stripe_id]);
% create set index
stripeIndex = unique([interp_struct.setID]);

stripeCounts = zeros(1,length(stripeIndex)); 

% initialize parpool
if exist('pool') ~= 1
    pool = parpool(pool_max);
end
% structure array to store the analysis data
outputs = struct;
local_meta = struct;
init_meta = struct;            

%Note: This was designed to allow for an array of AP groups to be passed
%in the same run request. On Savio it works better to create separate jobs
%for each bin. Keeping the structure as it shouldn't slow things down all
%that much and could be useful in the future
for K = state_vec
    for g = 1:length(stripe_groups)
        stripe_list = stripe_groups{g};
        % get list of all sets with traces conributing to current stripe
        stripe_sets = unique([interp_struct(ismember([interp_struct.stripe_id]...
            ,stripe_list)).setID]);
        for ss = 1:length(stripe_sets)
            removed_set = stripe_sets(ss);
            kept_sets = stripe_sets(stripe_sets~=removed_set);
            local_struct = struct;
            init_struct = struct;            
            mean_bin = mean(stripe_list);
            % initialize logL to - infinity
            logL_max = -Inf;

            % extract fluo_data
            trace_ind = find(ismember([interp_struct.stripe_sub_id],stripe_regions(1))...
                &ismember([interp_struct.stripe_id],stripe_list)...
                &ismember([interp_struct.setID],kept_sets));

            fluo_dat= cell([length(trace_ind), 1]);            
            for tr = 1:length(trace_ind)
                fluo_data{tr} = interp_struct(trace_ind(tr)).fluo;
            end            
            % random initialization of model parameters
            param_init = initialize_random (K, w, fluo_data);
            % approximate inference assuming iid data for param initialization
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, 1000, 1e-4);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);
            parfor i_local = 1:n_localEM
                % random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);

                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;
                noise_init = param_init.noise;

                init_struct(i_local).A_init = exp(A_log_init);                
                init_struct(i_local).v_init = v_init;
                init_struct(i_local).noise_init = noise_init;
%                 init_struct(i_local).set_id = s;
                init_struct(i_local).subset_id = i_local;
                % localEM call
                local_out = local_em_MS2_reduced_memory (fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);
                % Save Results 
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;
                local_struct(i_local).A_log = local_out.A_log;
                local_struct(i_local).A = exp(local_out.A_log);                
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;                
                lambda_inf = exp(local_out.lambda_log);
                local_struct(i_local).noise= 1/sqrt(lambda_inf);
                local_struct(i_local).pi0_log= local_out.pi0_log;
                local_struct(i_local).total_time = local_out.runtime;
                local_struct(i_local).total_steps = local_out.n_iter;
            end
            local_meta.init = init_struct;
            local_meta.local = local_struct;
            [logL, max_index] = max([local_struct.logL]);            
            outputs.local_runs = local_struct;
            outputs.pi0 =exp(local_struct(max_index).pi0_log);
            outputs.pi0_log = local_struct(max_index).pi0_log;

            outputs.v = local_struct(max_index).v(:);
            outputs.r = local_struct(max_index).r(:);

            outputs.noise = local_struct(max_index).noise;

            outputs.A = local_struct(max_index).A(:);
            outputs.A_mat = local_struct(max_index).A;
            outputs.A_log = local_struct(max_index).A_log;            
            
            outputs.stripe_id = stripe_groups{g};
            outputs.kept_sets = kept_sets;
            outputs.removed_set = removed_set;
            
            outputs.total_steps = local_struct(max_index).total_steps;
            outputs.total_time = local_struct(max_index).total_time;            
            
            outputs.set_size = length([interp_struct(trace_ind).time]);
            outputs.w = w;
            outputs.alpha = alpha;
            outputs.deltaT = deltaT;            
            
            fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) '_' date_str ...
                '_bin' num2str(stripe_list(1)) '_' num2str(stripe_list(end)) ...
                '_rm' num2str(removed_set) ];    
            
            out_file = [out_dir_single '/' fName_sub '_' num2str(now*100000) '.mat'];
            
            save(out_file, 'outputs');
        end           
    end
end
delete(pool);