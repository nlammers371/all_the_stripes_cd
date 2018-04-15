% Script to make aggregate activity profiles
addpath('../utilities/');
clear 
close all
%------------------------------Set System Params--------------------------%

%%% define ID variables
datatype = 'weka';
inference_type = 'set_bootstrap_results';
project = 'eve7stripes_inf_2018_03_27'; %project identifier
FigPath = ['../../fig/experimental_system/' project '/mRNA_maps/'];
mkdir(FigPath)
Tres = 20; % time resolution
%%% load data
load(['..\..\dat\' project '\inference_traces_' project '_dT20.mat']); % load inference traces
load(['..\..\dat\' project '\stripe_pos_' project '.mat']) % load stripe position info
load(['..\..\dat\' project '\fov_partitions_' project '.mat']);
%%% Define Mapping Variables
mRNA_decay_time = 7; %mRNA half life (minutes)
xDim_image = 1024; 
yDim_image = 256;
stripe_sep = 1024/64; % enforce fixed spacing between stripes mapping time
xDim_projection = 8*stripe_sep; % set scale of projection (intentionally sparse)
t_map = 35; % time used for mapping
plot_times = 1:50;
ap_start = .2*100;
ap_stop = .9*100;
[x_ref_mat, y_ref_mat] = meshgrid(1:xDim_image,1:yDim_image);
x_ref_image = 1:xDim_image; 
x_ref_projection = 1:xDim_projection;

%%% generate indexing vectors
stripe_index = [trace_struct_final.stripe_id_vec_interp];
stripe_index = unique(round(3*stripe_index(~isnan(stripe_index)))/3);
set_index = unique([trace_struct_final.setID]); % unique sets in data set
trace_time_vec = [];
trace_fluo_vec = [];
trace_x_vec = [];
trace_set_vec = [];
trace_stripe_vec = [];
trace_ap_vector = [];
% check for inverted sets
inversion_flags = zeros(1,length(set_index)); 
for i = 1:length(set_index)
    ap_vec = [trace_struct_final([trace_struct_final.setID]==i).ap_vector];
    x_vec = [trace_struct_final([trace_struct_final.setID]==i).xPos];
    x_ap_min = min(x_vec(ap_vec==min(ap_vec)));
    x_ap_max = max(x_vec(ap_vec==max(ap_vec)));
    inversion_flags(i) = x_ap_min>x_ap_max;
end

for i = 1:length(trace_struct_final)
    ft = trace_struct_final(i).fluo;
    tt = ceil(trace_struct_final(i).time/60);
    xp = trace_struct_final(i).xPos;
    setID = trace_struct_final(i).setID;
    if inversion_flags(set_index==setID)
        xp = xDim_image - xp + 1;
    end
    trace_stripe_vec = [trace_stripe_vec trace_struct_final(i).stripe_id_vec];
    trace_time_vec = [trace_time_vec tt(~isnan(ft))];
    trace_fluo_vec = [trace_fluo_vec ft(~isnan(ft))];
    trace_x_vec = [trace_x_vec xp];
    trace_set_vec = [trace_set_vec repelem(setID, length(trace_struct_final(i).xPos))];
    trace_ap_vector = [trace_ap_vector [trace_struct_final([trace_struct_final.setID]==i).ap_vector]];
end

%%% get average stripe positions. use to assign stripe positions
stripe_t_vec = stripe_pos_struct(1).t_vec;
map_filter = stripe_t_vec == t_map;
stripe_pos_mat= NaN(length(set_index),7);
for i = 1:length(set_index)
    stripe_id_mat = stripe_pos_struct(set_index(i)).stripe_id_mat(:,:,map_filter);
    if inversion_flags(i)
        stripe_id_mat = fliplr(stripe_id_mat);
    end
    set_stripes = unique(round(stripe_id_mat))';
    set_stripes = set_stripes(~isnan(set_stripes)); 
    set_stripes = set_stripes(set_stripes>0);
    ap_mat = fov_partitions(i).pixel_ap_id_mat;
    for s = 1:length(set_stripes)
        stripe_pos_mat(i,set_stripes(s)) = mean(ap_mat(stripe_id_mat==set_stripes(s)));
    end
end

stripe_positions_ap = nanmean(stripe_pos_mat); 
stripe_positions = (stripe_positions_ap-ap_start/100)*xDim_projection/(ap_stop-ap_start)*100;
xlabel_vec = ap_start:10:ap_stop;
xtick_vec = (xlabel_vec-ap_start)*xDim_projection/(ap_stop-ap_start);

%%% calculate mapping transformation for each set

% currently I apply a transformation to each set to generate a consistent
% mapping to the projection space. In future I will also account for
% "smearing" of x projection do to non-vertical stripes (eg 7)

set_map_filter_mat = NaN(length(set_index),xDim_image);
set_map_warp_factors = NaN(length(set_index),xDim_image);
for i = 1:length(set_index)
    stripe_id_mat = stripe_pos_struct(i).stripe_id_mat(:,:,map_filter);
    set_stripes = unique(round(stripe_id_mat))';
    set_stripes = set_stripes(~isnan(set_stripes)); 
    set_stripes = set_stripes(set_stripes>0);
    x_map_vector = NaN(1,xDim_image);    
    % mapping vectors
    x_orig_vec = []; 
    x_map_vec = [];
    % calculate x center for each stripe
    for j = set_stripes
        xp_stripe = round(mean(x_ref_mat(stripe_id_mat==j)));
        x_orig_vec = [x_orig_vec xp_stripe];
        x_map_vec = [x_map_vec stripe_positions(j)];
%         error('asfa')
    end    
    warp_factors = diff(x_map_vec)./diff(x_orig_vec); % mapping constants
    x_orig_centers = (x_orig_vec(2:end)+x_orig_vec(1:end-1))/2;
    offsets = (x_map_vec(2:end)+x_map_vec(1:end-1))/2 - ...
        x_orig_centers;
    offsets = offsets + x_orig_centers;
    x_warped_vec = (x_orig_vec(1:end-1)+x_orig_vec(2:end))/2.*warp_factors;
    
    set_map_filter_mat(i,x_ref_image<x_orig_vec(1)) = round(warp_factors(1).*(...
            x_ref_image(x_ref_image<x_orig_vec(1))-x_orig_centers(1))+offsets(1));
    set_map_warp_factors(i,x_ref_image<x_orig_vec(1)) = warp_factors(1);
    for j = 1:length(set_stripes)-1
        x_filter = x_ref_image<x_orig_vec(j+1)&x_ref_image>=x_orig_vec(j);        
        set_map_filter_mat(i,x_filter) = round(warp_factors(j).*(...
            x_ref_image(x_filter)-x_orig_centers(j))+offsets(j));
        set_map_warp_factors(i,x_filter) = warp_factors(j);
    end    
    set_map_filter_mat(i,x_ref_image>x_orig_vec(end)) = round(warp_factors(end).*(...
            x_ref_image(x_ref_image>x_orig_vec(end))-x_orig_centers(end)) + offsets(end));
    set_map_warp_factors(i,x_ref_image>x_orig_vec(end)) = warp_factors(end);         
end

%% ------------------ Make Instantaneous Production Maps ----------------%%
x_radius = 2; % radius of gauss kernel...nucleus diameter ~= 20-25
y_radius = 1;
kernel_sigma = 1; % this is kind of arbitrary
[x_ref_g, y_ref_g] = meshgrid(1:2*x_radius+1,1:2*y_radius+1);
x_ref_g = x_ref_g - x_radius - 1;
y_ref_g = y_ref_g - y_radius - 1;
r_mat = sqrt(x_ref_g.^2 + y_ref_g.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel

instant_mRNA_array = zeros(length(plot_times),xDim_projection,length(set_index)); % store profiles
stripe_pos_check_array = zeros(length(plot_times),xDim_projection,length(set_index)); % store profiles
for i = 1:length(set_index)
    set_filter = trace_set_vec == i;
    for t = plot_times
        time_filter = trace_time_vec==t;
        if sum(time_filter&set_filter) == 0
            continue
        end
        xp_vec = round(trace_x_vec(time_filter&set_filter));
        x_fluo_vec = trace_fluo_vec(time_filter&set_filter);
        x_stripe_vec = trace_stripe_vec(time_filter&set_filter);
        x_stripe_vec(x_stripe_vec~=round(x_stripe_vec)) = 0;
        x_map_vec = set_map_filter_mat(i,xp_vec); % mapping        
        xm_filter = x_map_vec>0;
        if sum(xm_filter) == 0
            continue
        end
        x_warp_vec = set_map_warp_factors(i,xp_vec); % warp adjustments
        x_sums = accumarray(x_map_vec(xm_filter)',(x_fluo_vec(xm_filter).*x_warp_vec(xm_filter))');
        x_sums_stripe = accumarray(x_map_vec(xm_filter)',x_stripe_vec(xm_filter)');
        x_index = unique(x_map_vec(xm_filter));
        x_sums = x_sums(ismember(1:max(x_map_vec),x_index));   
        x_sums_stripe = x_sums_stripe(ismember(1:max(x_map_vec),x_index));   
        instant_mRNA_array(t,x_index,i) = x_sums'; % save fluo_info        
        stripe_pos_check_array(t,x_index,i) = x_sums_stripe';                
    end
    norm_ref_array = ones(size(instant_mRNA_array,1),size(instant_mRNA_array,2));
%     norm_ref_array(instant_mRNA_array(:,:,i)==0) = 0;
    norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
    instant_mRNA_array(:,:,i) = conv2(instant_mRNA_array(:,:,i),g_kernel,'same');
    instant_mRNA_array(:,:,i) = instant_mRNA_array(:,:,i)./norm_ref_array;
    % set boundaries
    x_lower = min(set_map_filter_mat(i,round(trace_x_vec(set_filter))));
    x_upper = max(set_map_filter_mat(i,round(trace_x_vec(set_filter))));
    instant_mRNA_array(:,x_ref_projection<x_lower|x_ref_projection>x_upper,i) = NaN;
end
problem_sets = [];
set_vec = 1:11;
average_inst_mRNA = nanmean(instant_mRNA_array(:,:,~ismember(set_vec,problem_sets)),3);
average_stripe_pos = nanmean(stripe_pos_check_array,3);

%%% Make aggregate figure
instant_mRNA = figure;
colormap(jet(128));
imagesc(average_inst_mRNA)
h = colorbar;
ylabel(h,'mean fluorescence (AU)')
set(gca,'ytick',0:5:max(plot_times))
ylabel('minutes into nc14')
set(gca,'xtick',xtick_vec,'xticklabel',xlabel_vec)
xlabel('AP position')
title('Instantaneous Fluorescence Over Time (Aggregated)')
saveas(instant_mRNA,[FigPath 'instant_mRNA_all.png'],'png')
saveas(instant_mRNA,[FigPath 'instant_mRNA_all.pdf'],'pdf')

%%% make figs for single sets
for i = 1:length(set_index)
    instant_mRNA = figure;
    colormap(jet(128));
    imagesc(instant_mRNA_array(:,:,set_index(i)))
    h = colorbar;
    ylabel(h,'mean fluorescence (AU)')
    set(gca,'ytick',0:5:max(plot_times))
    ylabel('minutes into nc14')
    set(gca,'xtick',xtick_vec,'xticklabel',xlabel_vec)
    xlabel('AP position')
    title(['Instantaneous Fluorescence Over Time: Set ' num2str(set_index(i))])
    saveas(instant_mRNA,[FigPath 'instant_mRNA_' num2str(set_index(i)) '.png'],'png')
    saveas(instant_mRNA,[FigPath 'instant_mRNA_' num2str(set_index(i)) '.pdf'],'pdf')
end
close all
%%% ------------------ Make mRNA Accumulation Figure -------------------- %%
decay_kernel = fliplr(exp(-(plot_times-1)/mRNA_decay_time)); % convolution kernel to account for decay
accumulated_mRNA_array = zeros(length(plot_times),xDim_projection,length(set_index)); % store profiles
for i = 1:length(set_index)
    inst_mRNA = instant_mRNA_array(:,:,set_index(i)); % this is dumb
    for j = 1:size(inst_mRNA,2)
        for k = 1:size(inst_mRNA,1)            
            accumulated_mRNA_array(k,j,i) = sum(decay_kernel(end-k+1:end)'.*inst_mRNA(1:k,j));
        end
    end
end
average_acc_mRNA = nanmean(accumulated_mRNA_array,3);
%%% Make aggregate figure
accumulated_mRNA = figure;
colormap(jet(128));
imagesc(average_acc_mRNA)
h = colorbar;
ylabel(h,'accumualted mRNA (AU)')
set(gca,'ytick',0:5:max(plot_times))
ylabel('minutes into nc14')
set(gca,'xtick',xtick_vec,'xticklabel',xlabel_vec)
xlabel('AP position')
title('Accumulated mRNA Over Time (Aggregated)')
saveas(accumulated_mRNA,[FigPath 'accumualted_mRNA_all.png'],'png')
saveas(accumulated_mRNA,[FigPath 'accumulated_mRNA_all.pdf'],'pdf')
%%% make figs for single sets
for i = 1:length(set_index)
    accumulated_mRNA = figure;
    colormap(jet(128));
    imagesc(accumulated_mRNA_array(:,:,set_index(i)))
    h = colorbar;
    ylabel(h,'accumulated mRNA(AU)')
    set(gca,'ytick',0:5:max(plot_times))
    ylabel('minutes into nc14')
    set(gca,'xtick',xtick_vec,'xticklabel',xlabel_vec)
    xlabel('AP position')
    title(['Accumulated mRNA Over Time: Set ' num2str(set_index(i))])
    saveas(accumulated_mRNA,[FigPath 'accumulated_mRNA_' num2str(set_index(i)) '.png'],'png')
    saveas(accumulated_mRNA,[FigPath 'accumulated_mRNA_' num2str(set_index(i)) '.pdf'],'pdf')
end
close all
