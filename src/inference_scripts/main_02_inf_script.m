%Get environment variable from job script
stripe_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};
% stripe_groups = {1};
stripe_regions = [0];
%only do stripe centers for now
w = 8;
state_vec = [3,2];
%route to utilities folder
addpath('../utilities');
%------------------Define Inference Variables------------------------------%
% max num workers
pool_max = 24;
% set num local runs
n_localEM = 25;
% set max steps per inference
n_steps_max = 1000;
% set convergence criteria
eps = 10e-4;
datatype = 'weka';
project = 'eve7stripes_inf';

%Path to raw data
datapath = ['../../dat/' project '/'];
start_time = 0;
stop_time = 60;
% generate save names
dataname = ['inference_traces_w' num2str(w) '_' project '.mat'];
date_str = '2017-09-08_test';
out_dir =  ['../../out/' project '/' date_str '/' 'inference_w' ...
    num2str(w)  '/'];

if exist(out_dir) ~= 7
    mkdir(out_dir);
end
out_dir_single = [out_dir '/individual_results/']; 
if exist(out_dir_single) ~= 7
    mkdir(out_dir_single);
end
% if 1 inference will be conducted on "n_bootstrap" sets
bootstrap = 1;
n_bootstrap = 5;
sample_size = 6000;
% keep myself from doing stupid things....
if bootstrap == 0 && n_bootstrap ~= 1
    warning('Bootstrapping option not selected. Reseting n_bootstrap to 1')
    n_bootstrap = 1;
end

% Load data for inference into struct named: interp_struct
load([datapath dataname]);

%Get Tres and alpha
alpha = interp_struct(1).alpha;
deltaT = interp_struct(1).dT;

% create position index
stripeIndex = unique([interp_struct.stripe_id]);
stripeCounts = zeros(1,length(stripeIndex)); 
iter = 1;
for stripe = stripeIndex
    stripeCounts(iter) = sum([interp_struct([interp_struct.stripe_id]==stripe).N]);
    iter = iter + 1;
end

% initialize parpool
% pool = parpool(pool_max);
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
        for b = 1:n_bootstrap
            s = (g-1)*n_bootstrap + b;
            local_struct = struct;
            init_struct = struct;
            stripe_list = stripe_groups{g};
            mean_bin = mean(stripe_list);
            % initialize logL to - infinity
            logL_max = -Inf;

            % extract fluo_data
            trace_ind = find(ismember([interp_struct.stripe_sub_id],stripe_regions(1))...
                .*ismember([interp_struct.stripe_id],stripe_list));
            if bootstrap    
                set_size = sum([interp_struct(trace_ind).N]);
                ndp = 0;    
                sample_ids = [];
                if set_size < sample_size
                    warning('Set size less than prescribed inference sample size')
                end
                while ndp < sample_size
                    tr_id = randsample(trace_ind,1);
                    sample_ids = [sample_ids tr_id];
                    ndp = ndp + length(interp_struct(tr_id).time);
                end
                fluo_data = cell([length(sample_ids), 1]);
                for tr = 1:length(sample_ids)
                    fluo_data{tr} = interp_struct(sample_ids(tr)).fluo;
                end
            else
                fluo_data = cell([length(trace_ind), 1]);            
                for tr = 1:length(trace_ind)
                    fluo_data{tr} = interp_struct(trace_ind(tr)).fluo;
                end
            end
            % random initialization of model parameters
            param_init = initialize_random (K, w, fluo_data);
            % approximate inference assuming iid data for param initialization
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, 1000, 1e-4);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);
            for i_local = 1:n_localEM
                % random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);

                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;
                noise_init = param_init.noise;

                init_struct(i_local).A_init = exp(A_log_init);
                init_struct(i_local).R_init = logm(A_log_init)/deltaT;
                init_struct(i_local).v_init = v_init;
                init_struct(i_local).noise_init = noise_init;
                init_struct(i_local).set_id = s;
                init_struct(i_local).subset_id = i_local;
                % localEM call
                local_out = local_em_MS2_reduced_memory (fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);
                % Save Results 
                local_struct(i_local).set_id = s;
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;
                local_struct(i_local).A_log = local_out.A_log;
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).R = prob_to_rate(exp(local_out.A_log), deltaT);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / deltaT;
                lambda_inf = exp(local_out.lambda_log);
                local_struct(i_local).noise= 1/sqrt(lambda_inf);
                local_struct(i_local).pi0_log= local_out.pi0_log;
                local_struct(i_local).total_time = local_out.runtime;
                local_struct(i_local).total_steps = local_out.n_iter;
            end
            local_meta(s).init = init_struct;
            local_meta(s).local = local_struct;
            [logL, max_index] = max([local_struct.logL]);
            outputs(s).set_id = s;
            outputs(s).local_runs = local_struct;
            outputs(s).pi0 =exp(local_struct(max_index).pi0_log);
            outputs(s).pi0_log = local_struct(max_index).pi0_log;

            outputs(s).v = local_struct(max_index).v(:);
            outputs(s).r = local_struct(max_index).r(:);

            outputs(s).noise = local_struct(max_index).noise;

            outputs(s).A = local_struct(max_index).A(:);
            outputs(s).A_mat = local_struct(max_index).A;
            outputs(s).A_log = local_struct(max_index).A_log;

            outputs(s).R = local_struct(max_index).R(:);
            outputs(s).R_mat = local_struct(max_index).R;
            outputs(s).bin = mean_bin;
            outputs(s).bootstrap_flag = bootstrap;
            outputs(s).n_bootstrap = n_bootstrap;
            outputs(s).bootstrap_id = b;
            outputs(s).total_steps = local_struct(max_index).total_steps;
            outputs(s).total_time = local_struct(max_index).total_time;
            outputs(s).bootstrap_id = b;
            if bootstrap
                outputs(s).traces = sample_ids;                
                outputs(s).N = ndp;
            end
            outputs(s).set_size = length([interp_struct(trace_ind).time]);
            outputs(s).w = w;
            outputs(s).alpha = alpha;
            outputs(s).deltaT = deltaT;
            output = outputs(s);
            
            fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) '_' date_str '_bin' num2str(stripe_list(1)) '_' num2str(stripe_list(end)) '_group' num2str(b)];    
            suffix = 1;
            out_file = [out_dir_single '/' fName_sub   '_results_' num2str(suffix) '.mat'];
            while exist(out_file) == 2
                suffix = suffix + 1;
                out_file = [out_dir_single '/' fName_sub   '_results_' num2str(suffix) '.mat'];
            end
            save(out_file, 'output');
        end           
    end
    fName = ['eveSet_w' num2str(w) '_K' num2str(K) '_' date_str '_bin' num2str(stripe_list(1)) '_' num2str(stripe_list(end))];    
    suffix = 1;
    out_file = [out_dir '/' fName   '_results_' num2str(suffix) '.mat'];
    while exist(out_file) == 2
        suffix = suffix + 1;
        out_file = [out_dir '/' fName   '_results_' num2str(suffix) '.mat'];
    end
    save(out_file, 'outputs');
end
delete(pool)