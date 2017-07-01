addpath('../utilities');
% ---------------------- Generate synthetic data ----------------------
datapath = '../../processed_data/';
dataname = 'eveSet_2017_06_15.mat';
out_dir =  '../../inference_results';
outname = 'eveSet_2017_06_15';
if exist(out_dir) ~= 7
    mkdir(out_dir);
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
ap_groups = {[49,50], [51,52] , [53,54]};
K = 3;
alpha = interp_struct(1).alpha;
deltaT = interp_struct(1).dT;
w = 6 ; %interp_struct(1).w;
% max num workers
pool_max = 10;
% set num local runs
n_localEM = 25;
% set max steps per inference
n_steps_max = 1000;
% set convergence criteria
eps = 10e-4;
% initialize parpool
pool = parpool(12);
% structure array to store the analysis data
outputs = struct;
local_meta = struct;
init_meta = struct;

for s = 1:length(ap_groups)
    local_struct = struct;
    init_struct = struct;
    ap_list = ap_groups{s};
    mean_ap = mean(ap_list);
    % initialize logL to - infinity
    logL_max = -Inf;
    
    % extract fluo_data
    trace_ind = find(ismember([interp_struct.AP],ap_list));
    fluo_data = cell([length(trace_ind), 1]);
    t_pass = 1;
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
        init_struct(i_local).set_id = s;
        init_struct(i_local).subset_id = i_local;
        % localEM call
        local_out = local_em_MS2 (fluo_data, ...
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
        local_struct(i_local).total_time = local_out.total_time;
        local_struct(i_local).total_steps = local_out.n_iters;
    end
    local_meta(s).init = init_struct;
    local_meta(s).local = local_struct;
    [logL, max_index] = max([local_struct.logL]);
    outputs(s).set_id = s;
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
    outputs(s).AP = mean_ap;
    outputs(s).N = apCounts(ismember(apIndex,ap_list));
    outputs(s).w = w;
    outputs(s).alpha = alpha;
    outputs(s).deltaT = deltaT;
    outputs(s).total_time = local_struct(max_index).total_time;
    outputs(s).total_steps = local_struct(max_index).total_steps;
    output = outputs(s);
    local_meta_out = local_meta(s);

    fName = [outname '_ap' num2str(ap_list(1)) '_' num2str(ap_list(end))];
    % save the statistical validation results into a '.mat' file
    save([out_dir '/' fName '_results.mat'], 'output');
    %save([out_dir '/' fName '_all_inference_results.mat'], 'local_meta_out');
end

delete(pool)