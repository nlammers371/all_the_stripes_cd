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
nc_set_vec_ind = [nuclei_clean.setID];
tr_set_vec_ind = [trace_struct_final.setID];
set_index = unique(nc_set_vec_ind);
InterpGrid = trace_struct_final(1).InterpGrid;
peg_time = 49;
track_times = 25:round(peg_time/60);
stripe_pos_mat = NaN(length(track_times),7,length(set_index));
nc_pos_mat = NaN(length(track_times),7,length(set_index));
%%% define longform vectors to make scatters
nc_pos_vec_all = []; % longform vector to track each cohort nucleus
nc_set_vec_all = []; % corresponding set ID index
nc_time_vec_all = [];
nc_stripe_vec_all = [];
ncID_vec_all = [];

tr_pos_vec_all = []; % longform vector to track each cohort nucleus
tr_set_vec_all = []; % corresponding set ID index
tr_time_vec_all = [];
tr_stripe_vec_all = [];

[x_ref_mat,y_ref_mat] = meshgrid(1:1024,1:256); % position ref mats
% iterate through sets and save average stripe and nc cohort positions
stripe_plot_times = stripe_pos_struct(1).plot_times;
% stripe_plot_times = 25:50;%stripe_pos_struct(1).t_vec;

for i = 1:length(set_index)
    stripe_id_array = stripe_pos_struct(i).stripe_id_mat;
    set_nc = nuclei_clean(nc_set_vec_ind==set_index(i));
    set_tr = trace_struct_final(tr_set_vec_ind==set_index(i));
    % check for inverted sets
    tr_ap_vec = [set_tr.ap_vector_interp];
    tr_xp_vec = [set_tr.xPos_interp];
    tr_set_vec = repelem(set_index(i),length(tr_xp_vec));
    tr_time_vec = [set_tr.time_interp];
    tr_stripe_vec = [set_tr.stripe_id_vec_interp];
    
    nc_xp_vec = [set_nc.xPos_interp];
    nc_time_vec = [set_nc.time_interp];
    nc_set_vec = repelem(set_index(i),length(nc_time_vec));
    nc_stripe_vec = [set_nc.stripe_id_vec_interp];
    % this is really inefficient...    
    for j = 1:length(set_nc)
        ncID_vec = repelem(set_nc(j).ncID, length(set_nc(j).time_interp));
        ncID_vec_all = [ncID_vec_all ncID_vec];
    end
    
%     if tr_xp_vec(tr_ap_vec==min(tr_ap_vec)) > tr_xp_vec(tr_ap_vec==max(tr_ap_vec)) 
%         flipped = 1;
%         tr_xp_vec = 1024 - tr_xp_vec + 1;
%         nc_xp_vec = 1024 - nc_xp_vec + 1;
%     end        
    tr_pos_vec_all = [tr_pos_vec_all tr_xp_vec];
    tr_time_vec_all = [tr_time_vec_all tr_time_vec];
    tr_set_vec_all = [tr_set_vec_all tr_set_vec];
    tr_stripe_vec_all = [tr_stripe_vec_all tr_stripe_vec];
    
    nc_pos_vec_all = [nc_pos_vec_all nc_xp_vec];
    nc_time_vec_all = [nc_time_vec_all nc_time_vec];
    nc_set_vec_all = [nc_set_vec_all nc_set_vec];
    nc_stripe_vec_all = [nc_stripe_vec_all nc_stripe_vec];
end

% iterate through stripes and obtain cohort/stripe vectors
tr_stripe_struct = struct;
nc_stripe_struct = struct;
for s = 1:7
    % traces are simple
    tr_filter = tr_stripe_vec_all == s;
    tr_stripe_struct(s).tr_time_vec = tr_time_vec_all(tr_filter);
    tr_stripe_struct(s).tr_set_vec = tr_set_vec_all(tr_filter);
    tr_stripe_struct(s).tr_xp_vec = tr_pos_vec_all(tr_filter);
    % identify group of nuclei in stripe at appointed time
    nc_peg_filter = nc_stripe_vec_all==s & round(nc_time_vec_all/60)==peg_time;
    ncID_filter = ismember(ncID_vec_all,unique(ncID_vec_all(nc_peg_filter)));
    nc_stripe_struct(s).nc_time_vec = nc_time_vec_all(ncID_filter);
    nc_stripe_struct(s).nc_set_vec = nc_set_vec_all(ncID_filter);
    nc_stripe_struct(s).nc_xp_vec = nc_pos_vec_all(ncID_filter);
    nc_stripe_struct(s).nc_id_vec = ncID_vec_all(ncID_filter);
end
%% Make figures...
stripe_id = 1;
nc_t_vec = nc_stripe_struct(stripe_id).nc_time_vec/60;
nc_xp_vec = nc_stripe_struct(stripe_id).nc_xp_vec;
nc_setID_vec = nc_stripe_struct(stripe_id).nc_set_vec;
nc_setID_vec(nc_setID_vec~=8) = 1;
scatter(nc_xp_vec,nc_t_vec,45,nc_setID_vec)