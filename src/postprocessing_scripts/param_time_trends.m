% Compile results into summary structure. Generate summary plots
addpath('../utilities/');
% clear all
% close all
%------------------------------Set System Params--------------------------%
w = 7; %memory assumed for inference
K = 2; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
no_ends_flag = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 15;
%-----------------------------ID Variables--------------------------------%

% id variables
datatype = 'weka';
inference_type = 'set_bootstrap_results';
project = 'eve7stripes_inf_2018_02_20'; %project identifier
MarkerSize = 40; %Set size of markers for plots
plot_orig_rates = 0; %If 1 plot original rates (prior to fitting)
plot_scatters = 0; %If 1 plot single bootstrap results

%Generate filenames and writepath
% truncated_inference_w7_t20_alpha14_f1_cl1_no_ends1
id_string = [ 'truncated_inference_w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/t_window' num2str(t_window) '/' inference_type '/']; 

ReadPath = ['../../dat/' project '/' id_string];
FigPath = ['../../fig/experimental_system/' project '/' id_string];

mkdir(FigPath)

load([ReadPath 'hmm_results_mem' num2str(w) '_states' num2str(K) '.mat'])
t_vec = [hmm_results.t_inf];
inf_times = unique(t_vec);
region_vec = [hmm_results.binID];
inf_regions = unique(region_vec);

%%% Make parameter scatters
r_mat_mean = NaN(length(inf_times),length(inf_regions));
on_mat_mean = NaN(length(inf_times),length(inf_regions));
off_mat_mean = NaN(length(inf_times),length(inf_regions));

r_mat_ste = NaN(length(inf_times),length(inf_regions));
on_mat_ste = NaN(length(inf_times),length(inf_regions));
off_mat_ste = NaN(length(inf_times),length(inf_regions));

for i = 1:length(inf_times)
    for j = 1:length(inf_regions)
        r_mat_mean(i,j) = hmm_results(t_vec==inf_times(i)&region_vec==inf_regions(j)).initiation_mean(2)/60;
        r_mat_ste(i,j) = hmm_results(t_vec==inf_times(i)&region_vec==inf_regions(j)).initiation_std(2)/60;
        on_mat_mean(i,j) = [hmm_results(t_vec==inf_times(i)&region_vec==inf_regions(j)).R_fit_mean(2)];
        on_mat_ste(i,j) = [hmm_results(t_vec==inf_times(i)&region_vec==inf_regions(j)).R_fit_std(2)];
        off_mat_mean(i,j) = [hmm_results(t_vec==inf_times(i)&region_vec==inf_regions(j)).R_fit_mean(3)];
        off_mat_ste(i,j) = [hmm_results(t_vec==inf_times(i)&region_vec==inf_regions(j)).R_fit_std(3)];
    end
end
r_mat_ste = r_mat_ste / max(r_mat_mean(:));
r_mat_mean = r_mat_mean / max(r_mat_mean(:));
%% Make Scatter Figs   
ls = {};
for i = 1:7
    ls = [ls{:} {['stripe ' num2str(i)]}];
end
ls = [ls{:} {'full set'}];
s_path = [FigPath '/temp_scatters/'];
mkdir(s_path);
ref_vec = (2:22);
%%% r
inc = floor(128/7);
cm_indices = inc:inc:128;
cm = parula(128);

ref_cell = {ismember(ref_vec,2:3:22),ismember(ref_vec,3:3:22),ismember(ref_vec,4:3:22)};
title_cell = {'Anterior Stripe Flank', 'Stripe Center', 'Posterior Stripe Flank'};
save_cell = {'a','c','p'};
for i = 1:3
    r_scatter = figure;
    colormap(cm(cm_indices,:));
    hold on
%     h = colorbar;
%     set(h,'YTick',1:7)
    s = [];    
    for j = 1:7
        e = errorbar(inf_times/60,reshape(r_mat_mean(:,ref_vec==(3*j+i-2)),1,[]),...
            reshape(r_mat_ste(:,ref_vec==(3*j+i-2)),1,[]),'Color','black');
        e.CapSize = 0;
        s = [s scatter(inf_times/60,reshape(r_mat_mean(:,ref_vec==(3*j+i-2)),1,[]),...
            'MarkerFaceColor',cm((j-1)*inc+1,:), 'MarkerEdgeColor', 'black')];
    end
    e = errorbar(inf_times/60,reshape(r_mat_mean(:,1),1,[]),...
            reshape(r_mat_ste(:,1),1,[]),'Color','black');
    e.CapSize = 0;
    s = [s scatter(inf_times/60,reshape(r_mat_mean(:,1),1,[]),50,'s',...
            'MarkerFaceColor',[.3 .3 .3], 'MarkerEdgeColor', 'black')];
       
    legend(s,ls{:},'Location','southwest')
    xlabel('minutes into nc14')
    ylabel('normalized mRNA per second')
    axis([5 50 0 1.1*max(r_mat_mean(:))])
    title(['Estimated Initiation Rate (' title_cell{i} ')'])    
    saveas(r_scatter,[s_path '/r_scatter_' save_cell{i} '.png'],'png')
    saveas(r_scatter,[s_path '/r_scatter_' save_cell{i} '.pdf'],'pdf')
end


close all 
%%% kon
for i = 1:3
    k_on_scatter = figure;
    colormap(cm(cm_indices,:));
    hold on
%     h = colorbar;
%     set(h,'YTick',1:7)
    s = [];       
    for j = 1:7
        e = errorbar(inf_times/60,reshape(on_mat_mean(:,ref_vec==(3*j-(2-i))),1,[]),...
            reshape(on_mat_ste(:,ref_vec==(3*j-(2-i))),1,[]),'Color','black');
        e.CapSize = 0;
        s = [s scatter(inf_times/60,reshape(on_mat_mean(:,ref_vec==(3*j-(2-i))),1,[]),...
            'MarkerFaceColor',cm((j-1)*inc+1,:), 'MarkerEdgeColor', 'black')];        
    end
    e = errorbar(inf_times/60,reshape(on_mat_mean(:,1),1,[]),...
            reshape(on_mat_ste(:,1),1,[]),'Color','black');
    e.CapSize = 0;
    s = [s scatter(inf_times/60,reshape(on_mat_mean(:,1),1,[]),50,'s',...
            'MarkerFaceColor',[.3 .3 .3], 'MarkerEdgeColor', 'black')];
       
    legend(s,ls{:},'Location','northwest')
    xlabel('minutes into nc14')
    ylabel('switching rate (min^{-1})')
    axis([5 50 0 1.1*max(on_mat_mean(:))])
    title(['Estimated On Rate (' title_cell{i} ')'])    
    saveas(k_on_scatter,[s_path '/k_on_scatter_' save_cell{i} '.png'],'png')
    saveas(k_on_scatter,[s_path '/k_on_scatter_' save_cell{i} '.pdf'],'pdf')
end


%%% koff
for i = 1:3
    k_off_scatter = figure;
    colormap(cm(cm_indices,:));
    hold on
%     h = colorbar;
%     set(h,'YTick',1:7)
    s = [];
    for j = 1:7
        e = errorbar(inf_times/60,reshape(off_mat_mean(:,ref_vec==(3*j-(2-i))),1,[]),...
            reshape(off_mat_ste(:,ref_vec==(3*j-(2-i))),1,[]),'Color','black');
        e.CapSize = 0;
        s = [s scatter(inf_times/60,reshape(off_mat_mean(:,ref_vec==(3*j-(2-i))),1,[]),...
            'MarkerFaceColor',cm((j-1)*inc+1,:), 'MarkerEdgeColor', 'black')];
    end
    e = errorbar(inf_times/60,reshape(off_mat_mean(:,1),1,[]),...
            reshape(off_mat_ste(:,1),1,[]),'Color','black');
    e.CapSize = 0;
    s = [s scatter(inf_times/60,reshape(off_mat_mean(:,1),1,[]),50,'s',...
            'MarkerFaceColor',[.3 .3 .3], 'MarkerEdgeColor', 'black')];
    legend(s,ls{:},'Location','southwest')
    xlabel('minutes into nc14')
    ylabel('switching rate (min^{-1})')
    axis([5 50 0 1.1*max(off_mat_mean(:))])
    title(['Estimated Off Rate (' title_cell{i} ')'])    
    saveas(k_off_scatter,[s_path '/k_off_scatter_' save_cell{i} '.png'],'png')
    saveas(k_off_scatter,[s_path '/k_off_scatter_' save_cell{i} '.pdf'],'pdf')
end

%% Make Contour Maps
colormap(hot)
surf(on_mat_mean(1:end,:))
