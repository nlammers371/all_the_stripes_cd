% ---------------------- Generate synthetic data ----------------------
datapath = '../../processed_data/';
dataname = 'eveSet_2017_06_15.mat';
outpath = '../../inference_results';

if exist(outpath) ~= 7
    mkdir(outpath);
end
% Load data for inference into struct named: interp_struct
load([datapath dataname]);

% create AP index
apIndex = unique([interp_struct.AP]);
apCounts = zeros(1,length(apIndex)); 
iter = 1;
for ap = apIndex
    apCounts(iter) = sum([interp_struct([interp_struct.AP]==ap).N]);
    iter = iter + 1;
end
%------------------Define Inference Variables------------------------------%
ap_range = 40:41;
K = 3;
alpha = interp_struct(1).alpha;
deltaT = interp_struct(1).dT;
w = 2;
%w = interp_struct(1).w;
% initialize parpool
pool = parpool(10);
% structure to store synthetic data sets
data = struct;
% extract the synthetic data corresponding to the [set, i] pair
fluo_data = data(1).fluo_data;
% structure array to store the analysis data
outputs = struct;
local_meta = struct;
init_meta = struct;
i_iter = 1;

for set = 1:length(ap_range)
    local_struct = struct;
    init_struct = struct;
    ap = ap_range(set);
    % initialize logL to - infinity
    logL_max = -Inf;
    
    % extract fluo_data
    fluo_data = cell([length(trace_ind), 1]);
    tace_ind = find([interp_struct.AP]==ap);
    for tr = 1:length(trace_ind)
        fluo_data{tr} = interp_struct(trace_ind(tr)).fluo;
    end
    % random initialization of model parameters
    param_init = initialize_random (K, w, fluo_data);
    % approximate inference assuming iid data for param initialization
    local_iid_out = local_em_iid (fluo_data, param_init.v, ...
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
        init_struct(i_local).R_init = logm(A_log_init)/deltaT;
        init_struct(i_local).v_init = v_init;
        init_struct(i_local).noise_init = noise_init;
        init_struct(i_local).set_id = set;
        init_struct(i_local).subset_id = i_local;
        % localEM call
        local_out = local_em_MS2 (fluo_data, ...
            v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
            alpha, n_steps_max, eps);
        % Save Results 
        local_struct(i_local).set_id = set;
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
        local_struct(i_local).total_time = local_out.total_time;
        local_struct(i_local).total_steps = local_out.n_iters;
    end
    local_meta(set).init = init_struct;
    local_meta(set).local = local_struct;
    [logL, max_index] = max([local_struct.logL]);
    outputs(set).set_id = set;
    outputs(set).pi0 =exp(local_struct(max_index).pi0_log);
    outputs(set).pi0_log = local_struct(max_index).pi0_log;

    outputs(set).v = local_struct(max_index).v(:);
    outputs(set).r = local_struct(max_index).r(:);

    outputs(set).noise = local_struct(max_index).noise;

    outputs(set).A = local_struct(max_index).A(:);
    outputs(set).A_mat = local_struct(max_index).A;
    outputs(set).A_log = local_struct(max_index).A_log;

    outputs(set).R = local_struct(max_index).R(:);
    outputs(set).R_mat = local_struct(max_index).R;
    outputs(set).AP = ap;
    outputs(set).N = apCount(apIndex==ap);
    outputs(set).w = w;
    outputs(set).alpha = alpha;
    outputs(set).deltaT = deltaT;
    outputs(set).total_time = local_struct(max_index).total_time;
    outputs(set).total_steps = local_struct(max_index).total_steps;
end

% extract the current date in a string format
formatOut = 'yyyymmdd_HH_MM';
date_str = datestr(datetime('now'),formatOut);

% save the statistical validation results into a '.mat' file
save([out_dir '/' date_str '_results.mat'], 'outputs');
save([out_dir '/' date_str '_all_inference_results.mat'], 'local_meta');

% save the parameters used for data generation into a '.mat' file
save([dat_dir '/' date_str '_parameters.mat'], 'synthetic_parameters');

% save the generated data into a '.mat' file
save([dat_dir '/' date_str '_data.mat'], 'data');

% delete(pool)