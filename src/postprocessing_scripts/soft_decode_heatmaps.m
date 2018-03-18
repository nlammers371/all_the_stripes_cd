% Script to make figures illustrating spatio-temporal trends in hmm params
addpath('../utilities/');
% clear all
% close all
%%%%%%-----Set System Params
w = 7; %memory assumed for inference
K = 2; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = [bin_range_vec stripe_range(i) + j/3 - 2/3];
    end
end
bin_map_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_map_vec = [bin_map_vec stripe_range(i) + j*.3 - .6];
    end
end
% id variables
datatype = 'weka';
inference_type = 'set_bootstrap_results';
project = 'eve7stripes_inf_2018_02_20'; %project identifier

%Generate filenames and writepath
% truncated_inference_w7_t20_alpha14_f1_cl1_no_ends1
id_string = [ 'truncated_inference_w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins) '/' inference_type '/']; 

ReadPath = ['../../dat/' project '/' id_string];
FigPath = ['../../fig/experimental_system/' project '/' id_string '/surface_maps/'];
mkdir(FigPath)
% load results
load([ReadPath '/soft_decode_params.mat']);
t_start = 15; 
%---------------------- make surface maps -----------------------------%
%%
close all

time_vec = soft_decode_params.details(1).time_vec;
plot_times = t_start:5:max(time_vec)/60;
t_index_vec = 1:sum(time_vec/60>=t_start);

stripe_region_vec = soft_decode_params.stripe_region_vec;
s_index_vec = 1:length(stripe_region_vec);

% Initiation rate 
init_surf_fig = figure;
colormap(jet(128)/1.05)
init_rate_mat = soft_decode_params.initiation_rate_mat(time_vec/60>=t_start,:,2);
init_rate_mat = init_rate_mat / max(init_rate_mat(:));
imagesc(init_rate_mat, [0,1])
set(gca,'ytick',t_index_vec(ismember(time_vec(time_vec/60>=t_start)/60,plot_times)),'yticklabels',plot_times);
set(gca,'xtick',s_index_vec(ismember(stripe_region_vec,1:7)),'xticklabels',1:7)
title('Rate of Pol II Initiation')
h = colorbar;
ylabel(h, 'normalized mRNA per second')
saveas(init_surf_fig, [FigPath '/init_rate_surf.png'],'png')
saveas(init_surf_fig, [FigPath '/init_rate_surf.pdf'],'pdf')

% On rate 
on_surf_fig = figure;
colormap(jet(128)/1.05)
on_rate_mat = soft_decode_params.transition_prob_mat(time_vec/60>=t_start,:,2);
imagesc(on_rate_mat,[0 max(on_rate_mat(:))])
set(gca,'ytick',t_index_vec(ismember(time_vec(time_vec/60>=t_start)/60,plot_times)),'yticklabels',plot_times);
set(gca,'xtick',s_index_vec(ismember(stripe_region_vec,1:7)),'xticklabels',1:7)
title('Rate of Activation')
xlabel('eve stripe')
ylabel('time (minutes)')
h = colorbar;
ylabel(h, 'events per minute')
saveas(on_surf_fig, [FigPath '/on_rate_surf.png'],'png')
saveas(on_surf_fig, [FigPath '/on_rate_surf.pdf'],'pdf')

% Off rate 
off_surf_fig = figure;
colormap(jet(128)/1.05)
off_rate_mat = soft_decode_params.transition_prob_mat(time_vec/60>=t_start,:,3);
% off_rate_mat(off_rate_mat>2) = NaN;
imagesc(off_rate_mat,[0 max(off_rate_mat(:))])
set(gca,'ytick',t_index_vec(ismember(time_vec(time_vec/60>=t_start)/60,plot_times)),'yticklabels',plot_times);
set(gca,'xtick',s_index_vec(ismember(stripe_region_vec,1:7)),'xticklabels',1:7)
title('Rate of De-Activation')
xlabel('eve stripe')
ylabel('time (minutes)')
h = colorbar;
ylabel(h, 'events per minute')
saveas(off_surf_fig, [FigPath '/off_rate_surf.png'],'png')
saveas(off_surf_fig, [FigPath '/off_rate_surf.pdf'],'pdf')

% production rate 
p_surf_fig = figure;
colormap(jet(128)/1.05)
p_rate_mat = init_rate_mat.*(on_rate_mat./(on_rate_mat+off_rate_mat));
p_rate_mat = p_rate_mat / max(p_rate_mat(:));
imagesc(p_rate_mat, [0 1])
set(gca,'ytick',t_index_vec(ismember(time_vec(time_vec/60>=t_start)/60,plot_times)),'yticklabels',plot_times);
set(gca,'xtick',s_index_vec(ismember(stripe_region_vec,1:7)),'xticklabels',1:7)
title('Predicted Rate of Production')
xlabel('eve stripe')
ylabel('time (minutes)')
h = colorbar;
ylabel(h, 'normalized mRNA per second')
saveas(p_surf_fig, [FigPath '/predicted_production.png'],'png')
saveas(p_surf_fig, [FigPath '/predicted_production.pdf'],'pdf')

% production rate 
f_surf_fig = figure;
colormap(jet(128)/1.05)
f_rate_mat = soft_decode_params.mean_fluo_mat(time_vec/60>=t_start,:);
f_rate_mat = f_rate_mat / max(f_rate_mat(:));
% off_rate_mat(off_rate_mat>2) = NaN;
imagesc(f_rate_mat, [0 1])
set(gca,'ytick',t_index_vec(ismember(time_vec(time_vec/60>=t_start)/60,plot_times)),'yticklabels',plot_times);
set(gca,'xtick',s_index_vec(ismember(stripe_region_vec,1:7)),'xticklabels',1:7)
title('Actual Rate of Production')
xlabel('eve stripe')
ylabel('time (minutes)')
h = colorbar;
ylabel(h, 'normalized mRNA production rate')
saveas(f_surf_fig, [FigPath '/actual_production.png'],'png')
saveas(f_surf_fig, [FigPath '/actual_production.pdf'],'pdf')

p_diff_fig = figure;
colormap(jet(128)/1.05)
diff_mat = f_rate_mat - p_rate_mat;
imagesc(diff_mat,[-1,1]);
title('Actual Rate of Production')
xlabel('eve stripe')
ylabel('time (minutes)')
h = colorbar;
ylabel(h, 'difference btw predicted and actual mRNA production rates')
set(gca,'ytick',t_index_vec(ismember(time_vec(time_vec/60>=t_start)/60,plot_times)),'yticklabels',plot_times);
set(gca,'xtick',s_index_vec(ismember(stripe_region_vec,1:7)),'xticklabels',1:7)
saveas(p_diff_fig, [FigPath '/production_differential.png'],'png')
saveas(p_diff_fig, [FigPath '/production_differential.pdf'],'pdf')