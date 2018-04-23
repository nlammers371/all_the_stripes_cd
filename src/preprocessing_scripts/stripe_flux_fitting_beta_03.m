% Script to track stripe activity centers and (appx) mRNA levels across
% space and time
close all
clear 

% set filenames
project = 'eve7stripes_inf_2018_04_20'; %Project Identifier
fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory

flux_dynamics_path = [fig_path 'ap_positioning/stripe_dynamics/'];
mkdir(flux_dynamics_path);
trace_read_name = [data_path 'raw_traces_01_' project]; % names for compiled trace struct
trace_save_name = [data_path 'raw_traces_02_' project ]; % names for compiled trace struct
nucleus_read_name = [data_path 'ellipse_info_01_' project]; % names for compiled elipse struct
nucleus_save_name = [data_path 'ellipse_info_02_' project]; % names for compiled elipse struct
fov_name = [data_path 'fov_partitions.mat'];
cluster_name = [data_path 'stripe_clustering_results.mat'];
stripe_save_name = [data_path 'stripe_pos_' project '.mat'];
% load datasets
load(trace_read_name); % particle info
load(fov_name); % ap and stripe info at pixel level
load(nucleus_read_name);
load(cluster_name);

save_spline_figs = 1; % if 1 generates spline figs
%% generate fluorescence maps 
%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Color Info
cm = jet(128);
increment = floor(size(cm,1)/7);
stripe_colors = cm(1+((1:7)-1)*increment,:);
%%% get FOV info
xDim = size(fov_stripe_partitions(1).stripe_id_mat(:,:,1),2);
yDim = size(fov_stripe_partitions(1).stripe_id_mat(:,:,1),1);
%%% Smoothing Kernel
stripe_radius = round(mean([fov_stripe_partitions.ap_x_factor]))*1.5; % pixels 
kernel_radius = 60; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 15; 
[x_ref_kernel, y_ref_kernel] = meshgrid(1:2*kernel_radius+1,1:2*kernel_radius+1);
x_ref_kernel = x_ref_kernel - kernel_radius - 1;
y_ref_kernel = y_ref_kernel - kernel_radius - 1;
r_mat = sqrt(x_ref_kernel.^2 + y_ref_kernel.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel
g_kernel(r_mat>kernel_radius) = 0;
%%% fitting variables
t_window = 2; % lag/lead averaging window size in minutes
plot_times = fov_stripe_partitions(1).plot_times; % times during which to track stripe
min_time_all = 25; % first time to use for inference bins
min_time_5 = 30; % exception ro 5 as it comes on later
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% make spatial ref matrices
[x_ref_mat,y_ref_mat] = meshgrid(1:1024,1:256);
%%% make indexing vectors for traces
xPos_vec_particle = [trace_struct.xPos];
yPos_vec_particle = [trace_struct.yPos];
ap_vec_particle = [trace_struct.ap_vector];
set_vec_particle = [trace_struct.set_vector];
fluo_vec_particle = [trace_struct.fluo];
time_vec_particle = [trace_struct.time];
time_vec_particle = time_vec_particle(~isnan(fluo_vec_particle));
fluo_vec_particle = fluo_vec_particle(~isnan(fluo_vec_particle));
% same for nuclei
xPos_vec_nc = [schnitz_struct.xPos];
yPos_vec_nc = [schnitz_struct.yPos];
set_vec_nc = [schnitz_struct.set_vector];
time_vec_nc = [trace_struct.time];
set_index = unique(set_vec_particle); % indexing vector
stripe_class_vec = NaN(1,length(trace_struct));
stripe_pos_struct = struct; % save stripe location arrays
%%% perform stripe classifications
for i = 1:length(set_index)
    xp_set_vec = xPos_vec_particle(set_vec_particle==set_index(i));
    yp_set_vec = yPos_vec_particle(set_vec_particle==set_index(i));
    ap_set_vec = ap_vec_particle(set_vec_particle==set_index(i));
    fluo_set_vec = fluo_vec_particle(set_vec_particle==set_index(i));
    time_set_vec = time_vec_particle(set_vec_particle==set_index(i));    
    %%% extract stripe pixel map array
    stripe_mat_all = fov_stripe_partitions(i).stripe_id_mat;        
    %%% handle inversions    
    [~, ap_min] = min(ap_set_vec);
    [~, ap_max] = max(ap_set_vec);
    x_min = xp_set_vec(ap_min);
    x_max = xp_set_vec(ap_max);
    flip_flag = 0;
    if x_min > x_max                
        flip_flag = 1;
        xp_set_vec = size(stripe_mat_all,2) - xp_set_vec + 1;        
    end   
    stripe_id_vec_all = unique(stripe_mat_all(:))';
    if i == 3 % fix issue with s7 in 3
        stripe_id_vec_all(stripe_id_vec_all==7) = NaN;
    end
    stripe_id_vec_all = stripe_id_vec_all(~isnan(stripe_id_vec_all)); 
    temporal_fluo_array = NaN(size(stripe_mat_all));
    temporal_spline_mat_stable = NaN(size(stripe_mat_all,1),length(stripe_id_vec_all),length(plot_times));   
    temporal_spline_mat_raw = NaN(size(stripe_mat_all,1),length(stripe_id_vec_all),length(plot_times));   
    
    for m = 1:length(plot_times)
        j = length(plot_times) - m + 1;
        stripe_mat = stripe_mat_all(:,:,j);
        stripe_id_vec = unique(stripe_mat(:))';
        if i == 3 % don't fit stripe 7 in set 3
            stripe_id_vec(stripe_id_vec==7) = NaN;
        end
        stripe_id_vec = stripe_id_vec(~isnan(stripe_id_vec));
        stripe_frame_array = zeros(size(stripe_mat,1),size(stripe_mat,2));
        t = plot_times(j)*60;
        t_start = t - t_window*60;
        t_stop = t + t_window*60;
        t_fluo = fluo_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_x = xp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_y = yp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        idx = sub2ind(size(stripe_frame_array), t_y, t_x);
        stripe_frame_array(idx) = t_fluo;
        % apply smoothing kernel
        norm_ref_array = ones(size(stripe_frame_array));
        norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
        gauss_array = conv2(stripe_frame_array,g_kernel,'same');
        gauss_array = gauss_array./norm_ref_array;
        %%% matrices to store results        
        spline_mat_raw = NaN(size(stripe_mat,1),length(stripe_id_vec_all)); % spline fits
        spline_mat_stable = NaN(size(stripe_mat,1),length(stripe_id_vec_all)); % spline fits w/ back-fill
        for k = 1:length(stripe_id_vec)
            if stripe_id_vec(k) == 5 % account for fact that stripe 5 rises late
                min_time = min_time_5;
            else
                min_time = min_time_all;
            end
            stripe_id = stripe_id_vec(k);                                    
            wt_vec_raw = gauss_array(stripe_mat==stripe_id)';
            y_vec = y_ref_mat(stripe_mat==stripe_id)';
            indices = x_ref_mat(stripe_mat==stripe_id)';
            wt_vec = floor(wt_vec_raw/100);            
            if sum(wt_vec)==0
                continue
            end
            ind_wt_mat = repmat(indices,max(wt_vec),1);
            wt_mat = repmat(wt_vec,max(wt_vec),1);
            [~,y_grid] = meshgrid(1:size(wt_mat,2),1:size(wt_mat,1));
            ind_wt_mat(y_grid>wt_mat) = 0;                                        
            x_vec = reshape(ind_wt_mat,1,[]);
            y_vec = repelem(y_vec,max(wt_vec));
            y_vec = y_vec(x_vec>0);
            x_vec = x_vec(x_vec>0);           
            y_vec_sp = (1:size(stripe_mat,1))';
            pp = polyfit(y_vec,x_vec,3); % fit 4th degree polynomial to data
            poly = polyval(pp,y_vec_sp);             
            if  median(poly) <= stripe_radius || sum(abs(diff(poly))) > yDim...
                    ||median(poly) >= xDim - stripe_radius ||...
                    (i == 2 && plot_times(j)<35 && stripe_id_vec(k)==7) % special case for stripe 7 in set 2
                if plot_times(j) < min_time 
                       spline_mat_stable(:,ismember(stripe_id_vec_all,stripe_id_vec(k)))...
                            = nanmean(temporal_spline_mat_stable(:,ismember(stripe_id_vec_all,stripe_id_vec(k)),ismember(plot_times,min_time:min_time+4)),3);                        
                end
                continue
            end
            if plot_times(j) >= min_time
                spline_mat_stable(:,ismember(stripe_id_vec_all,stripe_id_vec(k))) = poly;                        
            else
                spline_mat_stable(:,ismember(stripe_id_vec_all,stripe_id_vec(k)))...
                    = nanmean(temporal_spline_mat_stable(:,ismember(stripe_id_vec_all,stripe_id_vec(k)),ismember(plot_times,min_time:min_time+4)),3);                        
            end
            spline_mat_raw(:,ismember(stripe_id_vec_all,stripe_id_vec(k))) = poly;               
        end                    
        if save_spline_figs
            poly_fig = figure;
            poly_fig.Visible = 'off';
            poly_fig.Position = [100 100 1024 256];
            hold on
            imagesc(gauss_array);           
            for k = 1:length(stripe_id_vec)
                stripe_id = stripe_id_vec(k);
                if stripe_id == 5
                    min_time = min_time_5;
                else
                    min_time = min_time_all;
                end
                if plot_times(j) <= min_time
                    plot(spline_mat_raw(:,ismember(stripe_id_vec_all,stripe_id_vec(k))),1:size(stripe_mat,1),'--','Color',cm(1+(stripe_id-1)*increment,:),...
                        'LineWidth',1)                
                    plot(spline_mat_stable(:,ismember(stripe_id_vec_all,stripe_id_vec(k))),1:size(stripe_mat,1),'Color',cm(1+(stripe_id-1)*increment,:),...
                        'LineWidth',2)
                else
                    plot(spline_mat_raw(:,ismember(stripe_id_vec_all,stripe_id_vec(k))),1:size(stripe_mat,1),'Color',cm(1+(stripe_id-1)*increment,:),...
                        'LineWidth',2)                
                end
            end
            axis([0 size(stripe_mat,2) 0 size(stripe_mat,1)])    
            ax = gca;
            ax.Visible = 'off';            
            mkdir([flux_dynamics_path '/set_' num2str(i) '/'])
            saveas(poly_fig,[flux_dynamics_path '/set_' num2str(i) '/ct_fluo_t' num2str(t) '.tif'],'tif');            
        end
        temporal_fluo_array(:,:,j) = stripe_frame_array;
        temporal_spline_mat_raw(:,:,j) = spline_mat_raw;        
        temporal_spline_mat_stable(:,:,j) = spline_mat_stable;        
    end        
    close all        
    %%%------------- Make time-dependent inference regions -------------%%%
    mean_center_mat = round(nanmean(temporal_spline_mat_stable(:,:,plot_times>min_time),3));
    stripe_id_mat_full = NaN(size(stripe_mat,1),size(stripe_mat,2),length(plot_times));         
    for t = 1:length(plot_times)        
        center_mat = round(temporal_spline_mat_stable(:,:,t));        
        neighbor_vec = [NaN nanmax(center_mat) NaN];
        active_indices = find(~isnan(nanmax(center_mat)));        
        for j = active_indices            
            for k = 1:size(stripe_mat,1)
                stripe_id_mat_full(k,max(1,center_mat(k,j)-stripe_radius):...
                                min(center_mat(k,j)+stripe_radius,xDim),t) = stripe_id_vec_all(j);
            end  
            if isnan(neighbor_vec(j))&& ~isnan(neighbor_vec(j+2)) % nan on left
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,max(1,center_mat(k,j) - 3*stripe_radius):...
                                max(center_mat(k,j)-stripe_radius-1,1),t) = stripe_id_vec_all(j) - 1/3;
                    stripe_id_mat_full(k,min(xDim,center_mat(k,j) + stripe_radius + 1):...
                                round(.5*(center_mat(k,j)+center_mat(k,j+1))),t) =...
                                stripe_id_vec_all(j) + 1/3;                            
                end                                                 
            elseif ~isnan(neighbor_vec(j))&& isnan(neighbor_vec(j+2)) % nan on right
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec_all(j) - 1/3;
                    stripe_id_mat_full(k,min(xDim,center_mat(k,j) + stripe_radius + 1):...
                                min(xDim,center_mat(k,j) + 3*stripe_radius),t) =...
                                stripe_id_vec_all(j) + 1/3;                             
                end     
            elseif isnan(neighbor_vec(j))&& isnan(neighbor_vec(j+2)) % nan on both sides
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,max(1,center_mat(k,j) - 3*stripe_radius):...
                                center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec_all(j) - 1/3;
                    stripe_id_mat_full(k,min(xDim,center_mat(k,j) + stripe_radius + 1):...
                                min(xDim,center_mat(k,j) + 3*stripe_radius),t) =...
                                stripe_id_vec_all(j) + 1/3;                             
                end     
            else
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                    center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec_all(j) - 1/3;
                    stripe_id_mat_full(k,min(xDim,center_mat(k,j) + stripe_radius + 1):...
                                    round(.5*(center_mat(k,j)+center_mat(k,j+1))),t) =...
                                    stripe_id_vec_all(j) + 1/3;                            
                end
            end
        end
    end  
    stripe_id_mat_show = stripe_id_mat_full;
    stripe_id_mat_show(isnan(stripe_id_mat_show)) = 0;
    mode_stripe_id_mat = mode(stripe_id_mat_show(:,:,plot_times>=min_time_5),3);    
    stripe_pos_struct(i).stripe_id_mat = stripe_id_mat_full;     
    stripe_RGB = NaN(size(stripe_mat,1),size(stripe_mat,2),3);    
    index_vec = 1:size(stripe_mat,2);        
    % make RGB image
    for j = 1:length(stripe_id_vec_all)
        stripe_id = stripe_id_vec_all(j);
        stripe_color = stripe_colors(stripe_id,:);
        for k = 1:3
            slice = stripe_RGB(:,:,k);
            slice(mode_stripe_id_mat==stripe_id-1/3) = stripe_color(k)/2;
            slice(mode_stripe_id_mat==stripe_id) = stripe_color(k)/1.5;
            slice(mode_stripe_id_mat==stripe_id+1/3) = stripe_color(k);
            slice(isnan(slice)) = .5;
            stripe_RGB(:,:,k) = slice;
        end
    end
    
    partition_fig = figure;  
    partition_fig.Visible = 'off';
    imshow(stripe_RGB)
    hold on
    for j = 1:length(stripe_id_vec_all)
        plot(mean_center_mat(:,j),1:size(stripe_mat,1),'Color','black','LineWidth',1.5)
    end
    title(['Set: ' num2str(i)])        
    saveas(partition_fig,[flux_dynamics_path '/inference_partitions_set' num2str(i) '.png'],'png')     
    %%% classify traces
    single_set_vec = [trace_struct.setID];
    set_indices = find(single_set_vec==i);
    
    for m = 1:length(set_indices)
        ind = set_indices(m);        
        xVec = round(trace_struct(ind).xPos);
        if flip_flag
            xVec = xDim - xVec + 1;
        end
        yVec = round(trace_struct(ind).yPos);
        t_trace = trace_struct(ind).time;
        f_trace = trace_struct(ind).fluo;
        t_trace = t_trace(~isnan(f_trace));
        tr_stripe_id_vec = zeros(1,length(t_trace))-1;
        fit_times = plot_times(plot_times>=min_time);
        for t = 1:length(t_trace)
            xp = xVec(t);
            yp = yVec(t);
            t_filter = round(t_trace(t)/60)==plot_times;
            if sum(t_filter) == 0 % skip early tp (will be back-filled)
                t_filter = 1; % set anything earlier than 10 to 10
            end
            s_id = stripe_id_mat_full(yp,xp,t_filter);            
            tr_stripe_id_vec(t) = s_id;                            
        end         
%         fill_ind = find(-1==(tr_stripe_id_vec));
%         id_vec = 1:length(tr_stripe_id_vec);        
%         t_ref_vec = plot_times;
%         t_ref_vec(t_ref_vec<min_time) = Inf;
%         first_ind = find(t_ref_vec>=min_time,1);
%         if ~isempty(fill_ind) 
%             % Assign early time points to nearest valid dynamic bin        
% %             for k = 1:length(fill_ind)
% %                 [~, nn] = min(abs(t_ref_vec - round(t_trace(fill_ind(k))/60))); % find nearest valid point
%                 tr_stripe_id_vec(fill_ind(k)) = stripe_id_mat_full(yVec(fill_ind(k)),xVec(fill_ind(k)),first_ind:first_ind+4);
%             end                
%         end        
        trace_struct(ind).stripe_id_inf = mode(tr_stripe_id_vec);
        trace_struct(ind).stripe_id_vec = tr_stripe_id_vec;
%         trace_struct(ind).xPos = xVec;
        trace_struct(ind).x_flip_flag = flip_flag;
        trace_struct(ind).stripe_add_flag  = 1;
    end    
    
    %%% classify nuclei
    single_set_vec = [schnitz_struct.setID];
    set_indices = find(single_set_vec==i);    
    for m = 1:length(set_indices)
        ind = set_indices(m);
        xVec = round(schnitz_struct(ind).xPos);
        yVec = round(schnitz_struct(ind).yPos);
        if flip_flag 
            xVec = xDim - xVec + 1;
        end
        t_trace = schnitz_struct(ind).time;                
        nc_stripe_id_vec = NaN(1,length(t_trace));
%         fit_times = plot_times(plot_times>=min_time);
        for t = 1:length(t_trace)
            xp = xVec(t);
            yp = yVec(t);
            t_filter = round(t_trace(t)/60)==plot_times;
            if sum(t_filter) == 0 %back-fill early tp
                t_filter = 1;
            end
            s_id = stripe_id_mat_full(yp,xp,t_filter);                        
            nc_stripe_id_vec(t) = s_id;                            
        end         
%         fill_ind = find(-1==(nc_stripe_id_vec));
%         id_vec = 1:length(nc_stripe_id_vec);        
%         t_ref_vec = plot_times;
%         t_ref_vec(t_ref_vec<min_time) = Inf;
%         if ~isempty(fill_ind) 
%             % Assign early time points to nearest valid dynamic bin        
%             for k = 1:length(fill_ind)
%                 [~, nn] = min(abs(t_ref_vec - round(t_trace(fill_ind(k))/60))); % find nearest valid point
%                 nc_stripe_id_vec(fill_ind(k)) = stripe_id_mat_full(yVec(fill_ind(k)),xVec(fill_ind(k)),nn);
%             end                
%         end        
        schnitz_struct(ind).stripe_id_inf = mode(nc_stripe_id_vec);
        schnitz_struct(ind).stripe_id_vec = nc_stripe_id_vec;
%         schnitz_struct(ind).xPos = xVec;
        schnitz_struct(ind).x_flip_flag = flip_flag;
        schnitz_struct(ind).stripe_add_flag = 1;
    end            
    stripe_pos_struct(i).plot_times = plot_times;
    disp(['Completed ' num2str(i) ' of ' num2str(length(set_index))])    
end
save(stripe_save_name, 'stripe_pos_struct')
save(trace_save_name,'trace_struct')
save(nucleus_save_name,'schnitz_struct')