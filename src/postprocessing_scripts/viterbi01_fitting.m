addpath('../utilities');
%memory assumed for inference
w = 7;
%states used for final inference
K = 3;
%set bin range
stripe_range = [1:7];
%write to csv?
write_csv = 0;
%Use aggregate fit rather than stripe-specific?
aggregate_fit = 1;
%Set time res
dT = 20;
start_time = 0;
stop_time = 60;
datatype = 'weka';
date_str = '2017-09-25';
project = 'eve7stripes_inf_2017_09_25';

% inpath = ['../../fig/experimental_system/' project '/' date_str '/mem' data_name '/hmm_results_3states_mem8.mat'];
datapath = ['../../fig/experimental_system\' project '\' date_str ...
    '\mem' num2str(w) '_states' num2str(K)  '\hmm_results_mem' num2str(w)...
    '_states' num2str(K) '.mat']; 
load(datapath)
if aggregate_fit
    OutPath = ['../../out/' project '/' date_str '/mem' ...
        num2str(w) '_states' num2str(K) '/'];
    ViterbiPath = ['../../out/' project '/' date_str '/mem' ...
        num2str(w) '_states' num2str(K) '/viterbi_fits_aggregate/' ];
else
    OutPath = ['../../out/' project '/' date_str '/mem' ...
        num2str(w) '_states' num2str(K) '/'];
    ViterbiPath = ['../../out/' project '/' date_str '/mem' ...
        num2str(w) '_states' num2str(K) '/viterbi_fits/' ];
end

if exist(ViterbiPath) ~= 7
    mkdir(ViterbiPath);
end

traces_all = load(['../../dat/' project '/inference_traces_t' num2str(dT) '_'...
    project '.mat']);
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.stripe_id]),stripe_range));
%Set header
header = {'SetID','ParticleID','StripeID','StripeSubID','Seconds', 'Fluo', ...
    'Fluo_Fit_Viterbi', 'Promoter_State_Viterbi','Stripe_Center', 'AP','x','y'};
%% Viterbi Fits
% w = traces_all(1).w;
%alpha (lenght of MS2 Loops in time steps)
alpha = hmm_results(1).alpha;
dT = hmm_results(1).dT;
% warning(ms2_loading_coeff (alpha, w))
%viterbi plots
job_size = 30;
s_iter = 1;
for bin = stripe_range
    if aggregate_fit
        hmm_bin = hmm_results([hmm_results.binID]==0);
    else
        hmm_bin = hmm_results([hmm_results.binID]==bin);
    end
    if isempty(hmm_bin)
        warning('No Inference Results for Stripe Region')
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
            v_fit.setID = bin_traces(ind).setID;
            v_fit.ParticleID = bin_traces(ind).OriginalParticle;
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
            save([ViterbiPath 'bin' num2str(bin) '_' 'trace_' num2str(ind)],'v_fit');
            % append viterbi and trace features to output array
            v_struct = bin_traces(ind);
            trace_stats = [double(v_struct.SetIDLong') double(v_struct.OrigParticleLong') double(v_struct.StripeIDLong') ...
                  double(v_struct.StripeSubIDLong') double(v_struct.time') double(v_struct.fluo') double(v_fit.fluo_viterbi') double(v_fit.z_viterbi')...
                  double(100*v_struct.StripeCenterLong') double(100*v_struct.ap_vector') double(v_struct.xPos') double(v_struct.yPos')];            
            if ind == 1 && s_iter == 1
                output_mat = trace_stats; 
            else
                output_mat = vertcat(output_mat, trace_stats); 
            end            
        end
    end
    s_iter = s_iter + 1;    
end

if aggregate_fit
    csvwrite_with_headers([OutPath '\' project '_dT' num2str(dT) '_longform_aggregate.csv'], ...
                       output_mat, header); 
else
    csvwrite_with_headers([OutPath '\' project '_dT' num2str(dT) '_longform.csv'], ...
                       output_mat, header); 
end               