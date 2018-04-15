% Script to Parse Relative Contributions from Bursting Parameters
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
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = [bin_range_vec stripe_range(i) + j/3 - 2/3];
    end
end

% id variables
datatype = 'weka';
inference_type = 'set';
project = 'eve7stripes_inf_2018_02_20'; %project identifier
FigPath = ['../../fig/experimental_system/' project '/mRNA_rate_profiles/'];
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
r_mat = [hmm_results.initiation_mean];
r_mat_ste = [hmm_results.initiation_std];
R_mat = [];
R_mat_ste = [];
for i = 1:length(hmm_results)
    R = hmm_results(i).R_fit_mean;
    R_ste = hmm_results(i).R_fit_mean;
    R_mat = [R_mat R(eye(K)~=1)'];
    R_mat_ste = [R_mat_ste R_ste(eye(K)~=1)'];
end
%%% generate data set for Mike

predicted_profile = r_mat(2,:).*(R_mat(1,:)./(R_mat(1,:) + R_mat(2,:)));
predicted_less_r = mean(r_mat(2,:))*(R_mat(1,:)./(R_mat(1,:) + R_mat(2,:)));
predicted_less_k_off = r_mat(2,:).*(R_mat(1,:)./(R_mat(1,:) + mean(R_mat(2,:))));
predicted_less_k_on = r_mat(2,:).*mean((R_mat(1,:))./(mean(R_mat(1,:)) + R_mat(2,:)));
predicted_only_k_on = mean(r_mat(2,:))*(R_mat(1,:)./(R_mat(1,:) + mean(R_mat(2,:))));

%%% generate empirical profile
fluo_profile = NaN(1,length(stripe_index));
t_start = t_inf-t_window/2;
t_stop = t_inf+t_window/2;
time_filter = trace_time_vec/60<t_stop & trace_time_vec/60>=t_start;
for i = 1:length(stripe_index)
    fluo_profile(i) = mean(trace_fluo_vec(time_filter&trace_stripe_vec==stripe_index(i)))/alpha_factor*60/Tres;
end

ss_param_mat = [double(stripe_index)',R_mat',R_mat_ste',r_mat',r_mat_ste',fluo_profile',...
    predicted_profile', predicted_less_r', predicted_less_k_on', predicted_less_k_off',...
    predicted_only_k_on'];

header = {'stripe_id', 'k_on_mean', 'k_off_mean','k_on_ste', 'k_off_ste', 'r_off_mean', 'r_on_mean',...
          'r_off_ste', 'r_on_ste','observed_production', 'predicted_production', 'predicted_less_r',...
          'predicted_less_k_on', 'predicted_less_k_off', 'predicted_only_k_on'};
csvwrite_with_headers([DataPath '\mean_param_profiles.csv'], ...
                       ss_param_mat, header,9); 

%%% make figures
cm = jet(128);
comparison_fig = figure;
hold on
a1 = area(stripe_index,fluo_profile,'FaceColor',cm(60,:));
a1.FaceAlpha = .4;
a2 = area(stripe_index,predicted_profile,'FaceColor',cm(20,:));
a2.FaceAlpha = .4;
xlim([min(stripe_index) max(stripe_index)]);
set(gca,'xtick',1:7)
xlabel('stripe')
ylabel('production rate (AU per min)')
title('Comparing Predicted and Actual Transcription Rates')
legend('actual','hmm prediction')
saveas(comparison_fig,[FigPath 'mRNA_rate_full.png'],'png')
saveas(comparison_fig,[FigPath 'mRNA_rate_full.pdf'],'pdf')

k_on_fig = figure;
hold on
a1 = area(stripe_index,fluo_profile,'FaceColor',cm(60,:));
a1.FaceAlpha = .4;
a2 = area(stripe_index,predicted_less_k_on,'FaceColor',cm(40,:));
a2.FaceAlpha = .4;
xlim([min(stripe_index) max(stripe_index)]);
set(gca,'xtick',1:7)
xlabel('stripe')
ylabel('production rate (AU per min)')
title('Effect of Removing  k_{on} Variation')
legend('actual','hmm prediction (less k_{on})')
saveas(k_on_fig,[FigPath 'mRNA_rate_less_k_on.png'],'png')
saveas(k_on_fig,[FigPath 'mRNA_rate_less_k_on.pdf'],'pdf')

k_on_fig2 = figure;
hold on
a1 = area(stripe_index,fluo_profile,'FaceColor',cm(60,:));
a1.FaceAlpha = .4;
a2 = area(stripe_index,predicted_only_k_on,'FaceColor',cm(40,:));
a2.FaceAlpha = .4;
xlim([min(stripe_index) max(stripe_index)]);
set(gca,'xtick',1:7)
xlabel('stripe')
ylabel('production rate (AU per min)')
title('Effect of k_{on} Alone (constant k_{off} and r)')
legend('actual','hmm prediction (only k_{on})')
saveas(k_on_fig2,[FigPath 'mRNA_rate_only_k_on.png'],'png')
saveas(k_on_fig2,[FigPath 'mRNA_rate_only_k_on.pdf'],'pdf')

k_off_fig = figure;
hold on
a1 = area(stripe_index,fluo_profile,'FaceColor',cm(60,:));
a1.FaceAlpha = .4;
a2 = area(stripe_index,predicted_less_k_off,'FaceColor',cm(80,:));
a2.FaceAlpha = .4;
xlim([min(stripe_index) max(stripe_index)]);
set(gca,'xtick',1:7)
xlabel('stripe')
ylabel('production rate (AU per min)')
title('Effect of Removing k_{off} Variation')
legend('actual','hmm prediction (less k_{off})')
saveas(k_off_fig,[FigPath 'mRNA_rate_less_k_off.png'],'png')
saveas(k_off_fig,[FigPath 'mRNA_rate_less_k_off.pdf'],'pdf')

r_fig = figure;
hold on
a1 = area(stripe_index,fluo_profile,'FaceColor',cm(60,:));
a1.FaceAlpha = .4;
a2 = area(stripe_index,predicted_less_r,'FaceColor',cm(110,:));
a2.FaceAlpha = .4;
xlim([min(stripe_index) max(stripe_index)]);
set(gca,'xtick',1:7)
xlabel('stripe')
ylabel('production rate (AU per min)')
title('Effect of Removing r Variation')
legend('actual','hmm prediction (less r)')
saveas(r_fig,[FigPath 'mRNA_rate_less_r.png'],'png')
saveas(r_fig,[FigPath 'mRNA_rate_less_r.pdf'],'pdf')