% Script to track stripe activity centers and (appx) mRNA levels across
% space and time
close all
clear all

% set filenames
project = 'eve7stripes_inf_2018_02_20'; %Project Identifier

fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory

flux_dynamics_path = [fig_path 'ap_positioning/stripe_dynamics/'];
mkdir(flux_dynamics_path);
trace_name = [data_path 'raw_traces_' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_' project]; % names for compiled elipse struct
fov_name = [data_path 'fov_partitions_' project '.mat'];
% load datasets
load(trace_name); % particle info
load(fov_name); % ap and stripe info at pixel level
save_spline_figs = 0; % if 1 generates spline figs
%% generate fluorescence maps 
cm = jet(128);
increment = floor(size(cm,1)/7);
stripe_radius = 35; % pixels .015;
kernel_radius = 60; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 35; % this is kind of arbitrary
[x_ref, y_ref] = meshgrid(1:2*kernel_radius+1,1:2*kernel_radius+1);
x_ref = x_ref - kernel_radius - 1;
y_ref = y_ref - kernel_radius - 1;
r_mat = sqrt(x_ref.^2 + y_ref.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel
% time stuff
t_window = 2; % lag/lead size in minutes
t_vec = 15:50;
t_start_partition = 15; % firts time to use for partitioning
% make indexing vectors
xPos_vec = [];
yPos_vec = [];
mean_xPos_vec = [];
mean_yPos_vec = [];
set_vec = [];
fluo_vec = [];
time_vec = [];
particle_vec = [];
for i = 1:length(trace_struct)
    fluo = trace_struct(i).fluo;
    time = trace_struct(i).time;
    time = time(~isnan(fluo));
    fluo = fluo(~isnan(fluo));    
    if isempty(fluo)
        error('asfa')
    end
    particle_vec = [particle_vec repelem(trace_struct(i).ParticleID,length(trace_struct(i).yPos))];
    fluo_vec = [fluo_vec fluo];
    time_vec = [time_vec time];
    xPos_vec = [xPos_vec trace_struct(i).xPos];
    yPos_vec = [yPos_vec trace_struct(i).yPos];
    mean_xPos_vec = [mean_xPos_vec mean(trace_struct(i).xPos(time>=600))];
    mean_yPos_vec = [mean_yPos_vec mean(trace_struct(i).yPos(time>=600))];
    set_vec = [set_vec repelem(trace_struct(i).setID, length(trace_struct(i).yPos))];    
end

set_index = unique(set_vec);
stripe_class_vec = NaN(1,length(trace_struct));

for i = 1:length(set_index)
    xp_set_vec = xPos_vec(set_vec==i);
    yp_set_vec = yPos_vec(set_vec==i);
    fluo_set_vec = fluo_vec(set_vec==i);
    time_set_vec = time_vec(set_vec==i);
    ap_mat = fov_partitions(i).pixel_ap_id_mat;
    stripe_mat = fov_partitions(i).pixel_stripe_id_mat;
    %%% handle inversions
    mean_ap_vec = nanmean(ap_mat);
    if mean_ap_vec(1) > mean_ap_vec(end)
        stripe_mat = fliplr(stripe_mat);
        ap_mat = fliplr(ap_mat);
        xp_set_vec = size(stripe_mat,2) - xp_set_vec + 1;
        mean_xPos_vec = size(stripe_mat,2) - mean_xPos_vec + 1;
    end
    stripe_id_vec = unique(reshape(stripe_mat(~isnan(stripe_mat)),1,[]),'stable');
    temp_fluo_array = NaN(size(stripe_mat,1),size(stripe_mat,2),length(t_vec));
    temp_spline_mat = NaN(size(stripe_mat,1),length(stripe_id_vec),length(t_vec));
    for j = 1:length(t_vec)
        frame_array = zeros(size(stripe_mat));
        t = t_vec(j)*60;
        t_start = t - t_window*60;
        t_stop = t + t_window*60;
        t_fluo = fluo_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_x = xp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_y = yp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        idx = sub2ind(size(frame_array), t_y, t_x);
        frame_array(idx) = t_fluo;
        norm_ref_array = ones(size(frame_array));
        norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
        gauss_array = conv2(frame_array,g_kernel,'same');
        gauss_array = gauss_array./norm_ref_array;
        centroid_mat = NaN(size(stripe_mat,1),length(stripe_id_vec)); % true centroids
        spline_mat = NaN(size(stripe_mat,1),length(stripe_id_vec)); % spline fits
        for k = 1:length(stripe_id_vec)
            stripe_id = stripe_id_vec(k);            
            for y = 1:size(stripe_mat,1)
                y_strip =  gauss_array(y,:);
                y_strip(stripe_mat(y,:)~=stripe_id) = 0;                
                indices = 1:size(stripe_mat,2);
                mi = sum(indices.*y_strip)/sum(y_strip);                
                centroid_mat(y,k) = mi;                
            end            
            center_vec = centroid_mat(:,k);
            y_vec = (1:size(stripe_mat,1))';
            f = fit(y_vec(~isnan(center_vec)),centroid_mat(~isnan(center_vec),k),'smoothingspline',...
                'SmoothingParam',0.01);            
            spline_mat(:,k) = feval(f,1:size(stripe_mat,1));
        end
        if save_spline_figs
            centroid_fig = figure;
            centroid_fig.Visible = 'off';
            centroid_fig.Position = [100 100 1024 256];
            hold on
            imagesc(gauss_array);
            scatter(reshape(centroid_mat,1,[]),repmat(1:size(stripe_mat,1),1,length(stripe_id_vec)),...
               5, 'MarkerFaceColor',[1 1 1]/4, 'MarkerEdgeColor',[1 1 1]/8)
            for k = 1:length(stripe_id_vec)
                stripe_id = stripe_id_vec(k);
                plot(spline_mat(:,k),1:size(stripe_mat,1),'Color',cm(1+(stripe_id-1)*increment,:),...
                        'LineWidth',1.5)
            end
            axis([0 size(stripe_mat,2) 0 size(stripe_mat,1)])
    %         text(10,25,[iIndex(round(t/60),2),...
    %         ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])
            ax = gca;
            ax.Visible = 'off';
            mkdir([flux_dynamics_path '/set_' num2str(i) '/'])
            saveas(centroid_fig,[flux_dynamics_path '/set_' num2str(i) '/ct_fluo_t' num2str(t) '.tif'],'tif');
        end
        temp_fluo_array(:,:,j) = frame_array;
        temp_spline_mat(:,:,j) = spline_mat;
    end    
    
    close all    
    %%% Make time-dependent inference regions
    mean_fluo_mat = nanmean(temp_fluo_array(:,:,t_vec>=t_start_partition),3);
    mean_center_mat = round(nanmean(temp_spline_mat(:,:,t_vec>=t_start_partition),3));
%     borders_array = NaN(size(mean_center_mat,1),size(mean_center_mat,2),4);
    stripe_id_mat_full = NaN(size(stripe_mat,1),size(stripe_mat,2),length(t_vec));     
    
    for t = 1:length(t_vec)        
        center_mat = round(temp_spline_mat(:,:,t));
        for j = 1:length(stripe_id_vec)
            for k = 1:size(stripe_mat,1)
                stripe_id_mat_full(k,center_mat(k,j)-stripe_radius:...
                                center_mat(k,j)+stripe_radius,t) = stripe_id_vec(j);
            end  
            if j == 1 %&& j ~= length(stripe_id_vec)
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,max(1,center_mat(k,j) - 3*stripe_radius):...
                                center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                round(.5*(center_mat(k,j)+center_mat(k,j+1))),t) =...
                                stripe_id_vec(j) + 1/3;                            
                end                                                 
            elseif j == length(stripe_id_vec)
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                min(1024,center_mat(k,j) + 3*stripe_radius),t) =...
                                stripe_id_vec(j) + 1/3;                             
                end     
            else
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                    center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                    round(.5*(center_mat(k,j)+center_mat(k,j+1))),t) =...
                                    stripe_id_vec(j) + 1/3;                            
                end
            end
        end
    end
    mode_stripe_id_mat = mode(stripe_id_mat_full(:,:,t_vec>=t_start_partition),3);
    mean_stripe_id_mat = mean(stripe_id_mat_full(:,:,t_vec>=t_start_partition),3);
    stripe_RGB = NaN(size(stripe_mat,1),size(stripe_mat,2),3);    
    index_vec = 1:size(stripe_mat,2);        
    % make RGB image
    for j = 1:length(stripe_id_vec)
        stripe_id = stripe_id_vec(j);
        stripe_color = cm(1+(stripe_id-1)*increment,:);
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
    for j = 1:length(stripe_id_vec)
        plot(mean_center_mat(:,j),1:size(stripe_mat,1),'Color','black','LineWidth',1.5)
    end
    title(['Set: ' num2str(i)])       
    saveas(partition_fig,[flux_dynamics_path '/inference_partitions_set' num2str(i) '.png'],'png') 
    
    %%% classify traces
    single_set_vec = [trace_struct.setID];
    set_indices = find(single_set_vec==i);
    
    for m = 1:length(set_indices)
        ind = set_indices(m);
        ParticleID = trace_struct(ind).ParticleID;
        xVec = xPos_vec(particle_vec==ParticleID);
        yVec = yPos_vec(particle_vec==ParticleID);
        t_trace = trace_struct(ind).time;
        f_trace = trace_struct(ind).fluo;
        t_trace = t_trace(~isnan(f_trace));
        stripe_id_vec = NaN(1,length(t_trace));
        for t = 1:length(t_trace)
            xp = xVec(t);
            yp = yVec(t);
            filter = round(t_trace(t)/60)==t_vec;
            if sum(filter) == 0
                continue
            end
            stripe_id_vec(t) = stripe_id_mat_full(yp,xp,filter);
        end         
        nan_ind = find(isnan(stripe_id_vec));
        if ~isempty(nan_ind) % at least 1 NaN            
            stripe_id_vec(nan_ind) = stripe_id_mat_full(round(mean(yVec(nan_ind))),round(mean(xVec(nan_ind))),1);
        end        
        trace_struct(ind).stripe_id_inf = mode(stripe_id_vec);
        trace_struct(ind).stripe_id_vec = stripe_id_vec;
    end    
    disp(['Completed ' num2str(i) ' of ' num2str(length(set_index))])
end

save(trace_name,'trace_struct')