% Script to call inference compilation function
%%% Define necessary variables
% Script to Conduct HMM Inference on Experimental Data
close all
clear 
addpath('E:\Nick\projects\hmmm\src\utilities'); % Route to hmmm utilities folder
% addpath('D:\Data\Nick\projects\hmmm\src\utilities'); % Route to hmmm utilities folder
savio = 0; % Specify whether inference is being conducted on Savio Cluster
ap_ref_index = 1:7;
ap_ref_index = reshape([ap_ref_index-1/3 ;ap_ref_index; ap_ref_index + 1/3],1,[]);

if savio
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    bin_groups = cell(1,length(savio_groups));
    for i =1:length(bin_groups)
        bin_groups{i} = ap_ref_index(savio_groups{i});
    end
else
    bin_groups = {};
    for i = 2:22
        bin_groups = [bin_groups{:} {round(i/3,1)}];
    end
end
%-------------------------------System Vars-------------------------------%
w = 7; % Memory
K = 2; % State(s) to use for inference
clipped_ends = 1; % if one, remove final w time steps from traces
dynamic_bins = 1; % if 1, use time-resolved region classifications
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_field = 1; % specify which fluo field to (1 or 3)
inference_times = 40*60;%(10:5:45)*60;
t_window = 30*60; % determines width of sliding window
inference_type = '_set';
project = 'eve7stripes_inf_2018_03_27_final';
Tres = 20;
DPFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
% call function
compile_inference_results(project,inference_type,w,K,Tres,alpha,fluo_type,clipped,...
    clipped_ends,dynamic_bins,t_window,t_inf,DPFolder)

