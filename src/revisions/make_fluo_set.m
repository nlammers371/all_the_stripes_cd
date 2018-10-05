% script to generate fluorescence-dependent binning for inference
clear 
close all
K = 3; % State(s) to use for inference
w = 7; % Memory
project = 'eve7stripes_inf_2018_04_28'; % project identifier
Tres = 20; % Time Resolution
% load trace data
datapath = ['../../dat/' project '/']; %Path to raw data
OutPath = '../../dat/revisions/';
mkdir(OutPath)
% generate read and write names
dataname = ['inference_traces_' project '_dT' num2str(Tres) '.mat'];
% Load data for inference into struct named: trace_struct_final
load([datapath dataname]);

n_bins = 10;
min_obs = 10;
% generate indexing vectors for grouping
mf_vec = NaN(1,numel(trace_struct_final));
stripe_vec = [trace_struct_final.stripe_id_inf];
pt_vec = [trace_struct_final.ParticleID];
for i = 1:numel(trace_struct_final)
    f_vec = trace_struct_final(i).fluo_interp;
    if numel(f_vec) >= min_obs
        mf_vec(i) = mean(f_vec);
    end
end

prctile_vec = linspace(100/n_bins, 100,n_bins);
f_id_vec_full = NaN(size(stripe_vec));
for i = 1:numel(prctile_vec)
    pct = prctile(mf_vec,prctile_vec(i));
    f_id_vec_full(mf_vec<=pct&isnan(f_id_vec_full)) = i;
end
i_pass = 1;
fluo_inf_struct = struct;
for i = 1:numel(pt_vec)    
    if ~isnan(mf_vec(i))
        fluo_inf_struct(i_pass).ParticleID = pt_vec(i);
        fluo_inf_struct(i_pass).stripe_id = stripe_vec(i);
        fluo_inf_struct(i_pass).FluoBin = f_id_vec_full(i);
        fluo_inf_struct(i_pass).fluo_interp = trace_struct_final(i).fluo_interp;
        fluo_inf_struct(i_pass).alpha_frac = trace_struct_final(i).alpha_frac;    
        i_pass = i_pass + 1;
    end
end
save([OutPath 'fluo_inf_struct.mat'] , 'fluo_inf_struct')