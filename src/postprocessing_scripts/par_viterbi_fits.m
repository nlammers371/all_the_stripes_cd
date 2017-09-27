addpath('../utilities');
%memory assumed for inference
w = 7;
%states used for final inference
K = 3;
%set bin range
bin_range = [-2:2];
dT = 20;
start_time = 0;
stop_time = 60;
datatype = 'weka';
date_str = '2017_09_25';
project = 'eve7stripes_inf_2017_09_26';
data_name = [ num2str(dT)];
% inpath = ['../../fig/experimental_system/' project '/' date_str '/mem' data_name '/hmm_results_3states_mem8.mat'];
datapath = ['../../fig/experimental_system\' project '\' date_str ...
    '\mem' num2str(w) '_states' num2str(K)  '\hmm_results_mem' num2str(w)...
    '_states' num2str(K) '.mat']; 
load(datapath)

OutPath = ['../../out/' project '/' date_str '/mem' ...
    num2str(w) '_states' num2str(K) '/viterbi_fits/' ];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

traces_all = load(['../../dat/' project '/inference_traces_t' data_name '_'...
    project '.mat']);
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.binID]),bin_range));
%initialize output array (will be written to csv)
output_mat = [];
header = {'SetID','ParticleID','StripeID','StripeSubID','Seconds', 'Fluo', ...
    'Fluo_Fit_Viterbi', 'Promoter_State_Viterbi','AP','x','y'};
%% Viterbi Plots
% w = traces_all(1).w;
%alpha (lenght of MS2 Loops in time steps)
alpha = hmm_results(1).alpha;
dT = hmm_results(1).dT;
% warning(ms2_loading_coeff (alpha, w))
%viterbi plots
job_size = 30;

for bin = stripe_range
    hmm_bin = hmm_results([hmm_results.binID]==bin);
    if isempty(hmm_bin)
        continue
    end
    v = hmm_bin.initiation_mean*dT;
    noise = hmm_bin.noise_mean;
    pi0_log = log(hmm_bin.pi0_mean);
    A_log = log(hmm_bin.A_mean);
    bin_traces = traces_all([traces_all.stripe_id]==bin);
    bin_traces = traces_all([bin_traces.stripe_sub_id]==0);
    viterbi_fits = struct;
    n_traces = length(bin_traces);
    n_jobs = ceil(n_traces/job_size);
    for n = 1:n_jobs
        viterbi_fits = struct;        
        batch_size = min(job_size,(n_traces-(n-1)*job_size));
        parfor i = 1:batch_size            
            ind = (n-1)*job_size+i;
            fluo = bin_traces(ind).fluo;            
            v_fit = viterbi (fluo, v, noise, pi0_log, ...
                                    A_log, K, w, alpha);
            v_fit.time_exp = bin_traces(ind).time;
            v_fit.fluo_exp = fluo;            
            v_fit.v = v;
            v_fit.w = w;
            v_fit.alpha = alpha;
            v_fit.noise = noise;
            v_fit.pi0 = exp(pi0_log);
            v_fit.A = exp(A_log);
            v_fit.trace_source = datapath;
            viterbi_fits(i).v_fit = v_fit;
    %         save([OutPath 'ap' num2str(ap) '_' 'trace_' num2str(i)],'v_fit');
        end
        for j = 1:length(viterbi_fits)
            v_fit = viterbi_fits(j).v_fit;
            ind = (n-1)*job_size + j;
            save([OutPath 'bin' num2str(bin) '_' 'trace_' num2str(ind)],'v_fit');
            % append viterbi and trace features to output array
            v_struct = traces_all(ind);
            trace_stats = [v_struct.SetIDLong' v_struct.OrigParticleLong' v_struct.StripeIDLong' ...
                  v_struct.StripeSubIDLong' v_struct.time' v_struct.fluo' v_fit.viterbi_fluo' v_fit.z_viterbi'...
                  v_struct.ap_vector' v_struct.xPos' v_struct.yPos'];
            output_mat = vertcat(output_mat, trace_stats);    
        end
    end
end
csvwrite_with_headers([inpath '\' dataname '_dT' num2str(dT) '_longform.csv'], ...
                       output_mat, header); 