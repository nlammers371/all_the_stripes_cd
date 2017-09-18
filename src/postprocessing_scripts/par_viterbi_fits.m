addpath('../utilities/');

%------------------------------Set System Params--------------------------%
%memory assumed for inference
w = 8;
%states used for final inference
K = 3;
%set bin range
stripe_range = [1:2];
sub_stripe_range = [0];
date_str = '2017-09-08_test';
project = ['eve7stripes_inf'];
% inpath = ['../../fig/experimental_system/' project '/' date_str '/mem' data_name '/hmm_results_3states_mem8.mat'];
inpath = ['D:\Data\Nick\projects\all_the_stripes_cd\fig\experimental_system\eve7stripes_inf\2017-09-08_test\mem8_states' num2str(K) '\'];
load([inpath '/hmm_results_' num2str(K) 'states_mem8.mat'])

OutPath = [inpath '/viterbi_fits/' ];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

traces_all = load('../../dat\eve7stripes_inf\inference_traces_w8_eve7stripes_inf.mat');
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.stripe_id]),stripe_range));

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
            v_fit.trace_source = '../../dat\mHMMeve2_weka_inf\inference_traces_w8_s0_35_mHMMeve2_weka_inf_old.mat';
            viterbi_fits(i).v_fit = v_fit;
    %         save([OutPath 'ap' num2str(ap) '_' 'trace_' num2str(i)],'v_fit');
        end
        for j = 1:length(viterbi_fits)
            v_fit = viterbi_fits(j).v_fit;
            num = (n-1)*job_size + j;
            save([OutPath 'bin' num2str(bin) '_' 'trace_' num2str(num)],'v_fit');
        end
    end
end