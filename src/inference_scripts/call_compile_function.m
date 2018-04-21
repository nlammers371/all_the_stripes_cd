% Script to call inference compilation function
%%% Define necessary variables
% Script to Conduct HMM Inference on Experimental Data
close all
clear 
% addpath('E:\Nick\projects\hmmm\src\utilities'); % Route to hmmm utilities folder
addpath('D:\Data\Nick\projects\hmmm\src\utilities'); % Route to hmmm utilities folder
ap_ref_index = 1:7;
ap_ref_index = reshape([ap_ref_index-1/3 ;ap_ref_index; ap_ref_index + 1/3],1,[]);


%-------------------------------System Vars-------------------------------%
w = 7; % Memory
K = 2; % State(s) to use for inference
clipped_ends = 1; % if one, remove final w time steps from traces
dynamic_bins = 1; % if 1, use time-resolved region classifications
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_type = 1; % specify which fluo field to (1 or 3)
inference_times = 25*60;%(10:5:45)*60;
t_window = 30; % determines width of sliding window
inference_type = 'dp';
project = 'eve7stripes_inf_2018_03_27_final';
Tres = 20;
alpha = 1.4;
t_inf = 40;
DPFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
% call function
compile_inference_results(project,inference_type,w,K,Tres,alpha,fluo_type,clipped,...
    clipped_ends,dynamic_bins,t_window,t_inf,DPFolder)


