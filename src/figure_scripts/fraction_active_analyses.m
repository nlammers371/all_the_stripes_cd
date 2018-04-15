% Script to estimate fraction of active nuclei as a function of space and
% time

addpath('../utilities/');
clear all
close all
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
project = 'eve7stripes_inf_2018_02_20'; %project identifier
OutPath = ['../../dat/' project '/' ];
FigPath = ['../../fig/experimental_system/' project '/fraction_active/'];
mkdir(OutPath)
mkdir(FigPath)
%---------------------------------Read in Files---------------------------%

%%% load inference traces 
load('D:\Data\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_02_20\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat');
load('D:\Data\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_02_20\inference_nuclei_eve7stripes_inf_2018_02_20_dT20.mat');
load('D:\Data\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_02_20\stripe_pos_eve7stripes_inf_2018_02_20.mat')

%%% Fraction on Analyses
trace_particle_vec = [trace_struct_final.ParticleID];
nc14_flag_vec = [nuclei_clean.nc14_flag];
analysis_nuclei = nuclei_clean(nc14_flag_vec==1);
nc_particle_vec = [analysis_nuclei.ParticleID];
% analysis_traces = trace_struct_final(ismember(trace_particle_vec,nc_particle_vec));
stripe_index = unique([trace_struct_final.stripe_id_inf]); 
stripe_index = stripe_index(~isnan(stripe_index));
InterpGrid = trace_struct_final(1).InterpGrid;
% make denominator (counts of nuclei per bin per time point)
nc_count_mat = zeros(length(InterpGrid),length(stripe_index));
tr_count_mat = zeros(length(InterpGrid),length(stripe_index));
for i = 1:length(analysis_nuclei)
    stripe_id_vec = analysis_nuclei(i).stripe_id_vec_interp;
    tt = analysis_nuclei(i).time_interp;
    ParticleID = analysis_nuclei(i).ParticleID;
    if ~isnan(ParticleID)
        tr_time = trace_struct_final(trace_particle_vec==ParticleID).time_interp;       
    end
    for j = 1:length(stripe_id_vec)
        nc_count_mat(InterpGrid==tt(j),stripe_index==stripe_id_vec(j)) = ...
            nc_count_mat(InterpGrid==tt(j),stripe_index==stripe_id_vec(j)) + 1;
        if ~isnan(ParticleID) 
            if sum(tr_time==tt(j)) == 1
                tr_count_mat(InterpGrid==tt(j),stripe_index==stripe_id_vec(j)) = ...
                tr_count_mat(InterpGrid==tt(j),stripe_index==stripe_id_vec(j)) + 1; 
            end
        end
    end             
end
%%
fraction_on_fig = figure;
colormap(jet(128))
imagesc(tr_count_mat./nc_count_mat)
colorbar;
title('Fraction Active Nuclei')
set(gca,'xtick',find(ismember(stripe_index,1:7)),'xticklabels',1:7)
set(gca,'ytick',find(ismember(InterpGrid,60*(0:5:50))),'yticklabels',0:5:50)
xlabel('stripe')
ylabel('minutes into nc14')
saveas(fraction_on_fig, [FigPath '/f_on_heatmap.png'],'png')
saveas(fraction_on_fig, [FigPath '/f_on_heatmap.pdf'],'pdf')
%% Find subset of traces with largest change in stripe ID
stripe_id_var_vec = [];
for i = 1:length(trace_struct_final)
    stripe_id_var_vec = [stripe_id_var_vec std(trace_struct_final(i).stripe_id_vec)];
end
%%
flux_traces = trace_struct_final(stripe_id_var_vec>.2);
for i = [27,28] %1:length(flux_traces)
    s_id_vec = flux_traces(i).stripe_id_vec;
    s_time = flux_traces(i).time;
    tr_time = flux_traces(i).time_interp;
    pID = flux_traces(i).ParticleID;
    fluo = flux_traces(i).fluo_interp;
    s_index = unique(s_id_vec);
    cm = jet(128);
    inc = floor(128/length(s_index));
    ls = {};
    s = [];
    blah_fig = figure;
    hold on
    for j = 1:length(s_index)
        st = min(s_time(s_id_vec==s_index(j)));
        ft = max(s_time(s_id_vec==s_index(j)));
        f_plot = fluo(tr_time>=st&tr_time<ft);
        t_plot = tr_time(tr_time>=st&tr_time<ft);
        plot(t_plot/60,f_plot/10e3,'Color','black')
        s = [s scatter(t_plot/60,f_plot/10e3,'MarkerFaceColor',cm(1 + (j-1)*inc,:),...
            'MarkerEdgeColor','black')];
        ls = [ls{:} {['s ' num2str(s_index(j))]}];
    end
    legend(s,ls{:})
    title(['Fluorescence as a Fucntion of Stripe Region (Particle ' num2str(pID)])
    xlabel('minutes into nc14')
    ylabel('au/10e3')
    tr_time = flux_traces(i).time_interp;
    saveas(blah_fig, [FigPath '/trace_flux_particle' num2str(pID) '.png'],'png')
    saveas(blah_fig, [FigPath '/trace_flux_particle' num2str(pID) '.pdf'],'pdf')
    
end
%% 
set_vec = [trace_struct_final.setID];
last_time_vec = zeros(1,11);
for i = 1:length(trace_struct_final)
    setID = trace_struct_final(i).setID;
    lt_curr = last_time_vec(setID);
    lt_new = max(lt_curr,max(trace_struct_final(i).time_interp));
    last_time_vec(setID) = lt_new;
end