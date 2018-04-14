% Script to examine divergence of stripe pattern and nuclei over time
close all
clear 
%%% set ID variables
project = 'eve7stripes_inf_2018_03_27_final';
DataPath = ['../../dat/' project];
FigPath = ['../../fig/' project '/Fig3/NucleusDynamics/'];
mkdir(FigPath)
Tres = 20; % t resolution used to generate data set
%%% load sets
load([DataPath '\inference_traces_' project '_dT' num2str(Tres) '.mat']); % load inference traces 
load([DataPath '\inference_nuclei_' project '_dT' num2str(Tres)  '.mat']); % load nuclei
load([DataPath '\stripe_pos_' project  '.mat']); % load inference traces 
%%% Define indexing vectors
nc_set_vec = [nuclei_clean.setID];
tr_set_vec = [trace_struct_final.setID];
set_index = unique(nc_set_vec);
InterpGrid = trace_struct_final(1).InterpGrid;
peg_time = 49*60;
track_times = 25:round(peg_time/60);
stripe_pos_mat = NaN(length(track_times),7,length(set_index));
nc_pos_mat = NaN(length(track_times),7,length(set_index));
[x_ref_mat,y_ref_mat] = meshgrid(1:1030,1:256); % position ref mats
% iterate through sets and save average stripe and nc cohort positions
stripe_plot_times = stripe_pos_struct(1).plot_times;
% stripe_plot_times = 25:50;%stripe_pos_struct(1).t_vec;

for i = 1:length(set_index)
    stripe_id_array = stripe_pos_struct(i).stripe_id_mat;
    set_nc = nuclei_clean(nc_set_vec==set_index(i));
    set_tr = trace_struct_final(tr_set_vec==set_index(i));
    % check for inverted sets
    ap_vec = [set_tr.ap_vector_interp];
    xp_vec = [set_tr.xPos_interp];
    flip_flag = 0;
    if xp_vec(ap_vec==min(ap_vec)) > xp_vec(ap_vec==max(ap_vec)) 
        flip_flag = 1;
    end
    cohort_vec = NaN(1,length(set_nc)); % keep track of nc identity
    for j = 1:length(set_nc)
        t_vec = round(set_nc(j).time_interp);
        stripe_vec = set_nc(j).stripe_id_vec_interp;
        t_filter = t_vec==peg_time;
        if sum(t_filter) ~= 0
            cohort_vec(j) = stripe_vec(t_filter);
        end
    end
    stripe_index = unique(round(cohort_vec));
    stripe_index = stripe_index(~isnan(stripe_index));    
    ref_pos_nc = NaN(size(stripe_index));
    ref_pos_stripe = NaN(size(stripe_index));
    for k = 1:length(track_times)
        t = length(track_times) - k + 1;
        tt = track_times(t);
        for s = 1:length(stripe_index)
            mean_x_stripe = mean(x_ref_mat(stripe_id_array(:,:,tt==stripe_plot_times)==stripe_index(s)));
            if k == 1
                ref_pos_stripe(s) = mean_x_stripe;
            end            
            stripe_pos_mat(t,stripe_index(s),i) = mean_x_stripe - ref_pos_stripe(s);
            cohort_list = find(cohort_vec == stripe_index(s));
            xp_vec_nc = [];
            for c = cohort_list
                t_vec = round(set_nc(c).time_interp/60);
                if sum(t_vec==tt) > 0
                    if flip_flag == 1
                        xp = 1024 - nanmean(set_nc(c).xPos_interp(t_vec==tt)) + 1;
                    else
                        xp = nanmean(set_nc(c).xPos_interp(t_vec==tt));
                    end
                    xp_vec_nc = [xp_vec_nc xp];
                end
            end
            if k == 1
                ref_pos_nc(s) = mean(xp_vec_nc);
            end
            nc_pos_mat(t,stripe_index(s),i) = mean(xp_vec_nc) - ref_pos_nc(s);
        end
    end    
end
%% Make Combined Fig
cm = jet(128);
increment = floor(size(cm,1)/7);
stripe_colors = cm(1+((1:7)-1)*increment,:);

nc_pos_mean = nanmean(nc_pos_mat,3);
stripe_pos_mean = nanmean(stripe_pos_mat,3);

cohort_fig = figure;
hold on
for i = 1:7
    plot(track_times,nc_pos_mean(:,i),'--','Color',stripe_colors(i,:),'LineWidth',1.5)
    plot(track_times,stripe_pos_mean(:,i),'Color',stripe_colors(i,:),'LineWidth',1.5)
end
plot(track_times,repelem(0,length(track_times)),'Color','black')
xlabel('minutes into nc14')
ylabel('pixels')