% Script to make plots of time-averaged HMM parameters
addpath('../utilities/');
clear 
close all
%%%%%%-----Set System Params
w = 7; %memory assumed for inference
K = 2; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
fluo_field = 1;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 30;
t_inf = 35;
t_start_stripe = 10; 
%-----------------------------ID Variables--------------------------------%
% id variables
datatype = 'weka';
inference_type = 'set';
project = 'eve7stripes_inf_2018_02_20'; %project identifier
FigPath = ['../../fig/experimental_system/' project '/ss_hmm_parameters/'];
mkdir(FigPath)
%Generate filenames and writepath
id_thing = [ '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins)  '/']; 

DataPath = ['../../dat/' project '/' id_thing '/K' num2str(K) '_summary_stats/' ];
mkdir(DataPath)
% load inference summary 
load([DataPath 'hmm_results_mem' num2str(w) '_states' num2str(K) '_ti' num2str(t_inf)...
    '_tw' num2str(t_window) '.mat'])
% load traces
load(['..\..\dat\' project '\inference_traces_' project '_dT' num2str(Tres) '.mat']);

%%% generate reference vectors 
alpha_factor = w-alpha*.5; % account for MS2 rise time
trace_stripe_vec = round(3*[trace_struct_final.stripe_id_vec_interp])/3;
stripe_index = unique(trace_stripe_vec);
stripe_index = stripe_index(~isnan(stripe_index));
trace_fluo_vec = [trace_struct_final.fluo_interp];
trace_time_vec = [trace_struct_final.time_interp];
%%% generate comparison profiles

%%% aggregate inf results
r_mat_mean = [hmm_results.initiation_mean];
r_mat_ste = [hmm_results.initiation_std];

R_mat_mean = [];
R_mat_ste = [];
for i = 1:length(hmm_results)
    R = hmm_results(i).R_fit_mean;
    R_ste = hmm_results(i).R_fit_std;
    R_mat_mean = [R_mat_mean R(eye(K)~=1)'];
    R_mat_ste = [R_mat_ste R_ste(eye(K)~=1)'];
end
R_mat_mean = flipud(R_mat_mean);
R_mat_ste = flipud(R_mat_ste);
cm = jet(128);
%%% Make figures
r_ss_fig = figure;
hold on
for k = 1:K
    e = errorbar(stripe_index,r_mat_mean(k,:),r_mat_ste(k,:),'Color','black');
    e.CapSize = 0;
    scatter(stripe_index,r_mat_mean(k,:),'MarkerFaceColor',cm(30*k,:),'MarkerEdgeColor','black')
end
title('Time-Averaged Initiation Rates')
xlabel('stripe')
set(gca,'xTick',1:7)
ylabel('AU per minute')
saveas(r_ss_fig,[FigPath 'ss_initiation_rates.png'],'png')
saveas(r_ss_fig,[FigPath 'ss_initiation_rates.pdf'],'pdf')

rate_ss_fig = figure;
hold on
for k = (1:K)
    e = errorbar(stripe_index,R_mat_mean(k,:),R_mat_ste(k,:),'Color','black');
    e.CapSize = 0;
    scatter(stripe_index,R_mat_mean(k,:),'MarkerFaceColor',cm(30*k,:),'MarkerEdgeColor','black')
end
title('Time-Averaged Transition Rates')
xlabel('stripe')
set(gca,'xTick',1:7)
ylabel('events per minute')
saveas(rate_ss_fig,[FigPath 'ss_transition_rates.png'],'png')
saveas(rate_ss_fig,[FigPath 'ss_transition_rates.pdf'],'pdf')