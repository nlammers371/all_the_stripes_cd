% Script to track stripe activity centers and (appx) mRNA levels across
% space and time
close all
clear 

% set filenames
project = 'eve7stripes_inf_2018_03_27_beta'; %Project Identifier

fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory

flux_dynamics_path = [fig_path 'ap_positioning/stripe_dynamics/'];
mkdir(flux_dynamics_path);
trace_name = [data_path 'raw_traces_' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_' project]; % names for compiled elipse struct
fov_name = [data_path 'fov_partitions_' project '.mat'];
stripe_pos_name = [data_path 'stripe_boundaries_' project '.mat'];
stripe_save_name = [data_path 'stripe_pos_' project '.mat'];
% load datasets
load(trace_name); % particle info
load(fov_name); % ap and stripe info at pixel level
load(nucleus_name);
load(stripe_pos_name);
save_stripe_figs = 1; % if 1 generates spline figs
%% track stripe dynamics over time
%%% Color Info
cm = jet(128);
increment = floor(size(cm,1)/7);

%%% fit variables
max_disp = 60; % max permissible fluo dispersion for stripe classfication 
stripe_radius = 35; % pixels .015;
kernel_radius = 30; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 15; % this is kind of arbitrary
[x_ref, y_ref] = meshgrid(1:2*kernel_radius+1,1:2*kernel_radius+1);
x_ref = x_ref - kernel_radius - 1;
y_ref = y_ref - kernel_radius - 1;
r_mat = sqrt(x_ref.^2 + y_ref.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel
g_kernel(r_mat>kernel_radius) = 0; % make circular

%%% time stuff
t_critical = fov_partitions(1).center_time; % get center time...this determines starting point for inference
t_window = 2; % lag/lead averaging window size in minutes
% t_vec = 25:50; %start fitting stripes at 25 min
t_fit_vec = 1:50;
%%% make indexing vectors for traces
xPos_vec_particle = [];
yPos_vec_particle = [];
mean_xPos_vec_particle = [];
mean_yPos_vec_particle = [];
set_vec_particle = [];
fluo_vec_particle = [];
time_vec_particle = [];
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
    fluo_vec_particle = [fluo_vec_particle fluo];
    time_vec_particle = [time_vec_particle time];
    xPos_vec_particle = [xPos_vec_particle trace_struct(i).xPos];
    yPos_vec_particle = [yPos_vec_particle trace_struct(i).yPos];
    mean_xPos_vec_particle = [mean_xPos_vec_particle mean(trace_struct(i).xPos(time>=600))];
    mean_yPos_vec_particle = [mean_yPos_vec_particle mean(trace_struct(i).yPos(time>=600))];
    set_vec_particle = [set_vec_particle repelem(trace_struct(i).setID, length(trace_struct(i).yPos))];    
end
% same for nuclei
xPos_vec_nc = [];
yPos_vec_nc = [];
mean_xPos_vec_nc = [];
mean_yPos_vec_nc = [];
set_vec_nc = [];
fluo_vec_nc = [];
time_vec_nc = [];
nc_vec = [];
for i = 1:length(schnitz_struct)    
    time = schnitz_struct(i).time;                
    nc_vec = [nc_vec repelem(schnitz_struct(i).ncID,length(schnitz_struct(i).yPos))];
    time_vec_nc = [time_vec_nc time];
    xPos_vec_nc = [xPos_vec_nc schnitz_struct(i).xPos];
    yPos_vec_nc = [yPos_vec_nc schnitz_struct(i).yPos];
    mean_xPos_vec_nc= [mean_xPos_vec_nc mean(schnitz_struct(i).xPos(time>=600))];
    mean_yPos_vec_nc = [mean_yPos_vec_nc mean(schnitz_struct(i).yPos(time>=600))];
    set_vec_nc = [set_vec_nc repelem(schnitz_struct(i).setID, length(schnitz_struct(i).yPos))];    
end

set_index = unique(set_vec_particle);
stripe_class_vec = NaN(1,length(trace_struct));
stripe_pos_struct = struct; % save stripe location arrays

%%
for i = 1:length(set_index)
    xp_set_vec = xPos_vec_particle(set_vec_particle==set_index(i));
    yp_set_vec = yPos_vec_particle(set_vec_particle==set_index(i));
    fluo_set_vec = fluo_vec_particle(set_vec_particle==set_index(i));
    time_set_vec = time_vec_particle(set_vec_particle==set_index(i));
    ap_mat = fov_partitions(i).pixel_ap_id_mat;
    stripe_mat = fov_partitions(i).pixel_stripe_id_mat;
    %%% handle inversions
    mean_ap_vec = nanmean(ap_mat);
    if mean_ap_vec(1) > mean_ap_vec(end)
        stripe_mat = fliplr(stripe_mat);
        ap_mat = fliplr(ap_mat);
        xp_set_vec = size(stripe_mat,2) - xp_set_vec + 1;
        mean_xPos_vec_particle = size(stripe_mat,2) - mean_xPos_vec_particle + 1;
    end
    stripe_id_vec = unique(reshape(stripe_mat(~isnan(stripe_mat)),1,[]),'stable');
    temp_fluo_array = NaN(size(stripe_mat,1),size(stripe_mat,2),length(t_fit_vec));
    temp_center_mat = NaN(size(stripe_mat,1),length(stripe_id_vec),length(t_fit_vec));
    temp_disp_mat = NaN(size(stripe_mat,1),length(stripe_id_vec),length(t_fit_vec));
    
    %%% first fit stripe position at critical time
    t_ind = find(t_fit_vec==t_critical/60);
    tp = t_critical; 
    x_tol = 1.5*stripe_radius; % allow flexibility for this piece
    initial_centers = stripe_positions{i};
    frame_array = zeros(size(stripe_mat)); 

    t_start = tp - t_window*60;
    t_stop = tp + t_window*60;

    t_fluo = fluo_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
    t_x = xp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
    t_y = yp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
    idx = sub2ind(size(frame_array), t_y, t_x);
    frame_array(idx) = t_fluo;
    % apply smoothing kernel
    norm_ref_array = ones(size(frame_array));
    norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
    gauss_array = conv2(frame_array,g_kernel,'same');
    gauss_array = gauss_array./norm_ref_array;    
    
    %%% loop through stripes in set
    for k = 1:length(stripe_id_vec)
        stripe_id = stripe_id_vec(k);                          
        st_gauss =  gauss_array;
        st_gauss(stripe_mat~=stripe_id) = 0;   
        st_gauss_vec = reshape(st_gauss,1,[]);
        x_indices = repelem(1:size(stripe_mat,2),size(stripe_mat,1));                    
        y_indices = repmat(1:size(stripe_mat,1),1,size(stripe_mat,2));                    
        %%% apply weights
        wt_vec = floor(st_gauss_vec/100);
        if max(wt_vec) == 0
            continue
        end
        x_ind_wt_mat = repmat(x_indices,max(wt_vec),1);
        y_ind_wt_mat = repmat(y_indices,max(wt_vec),1);
        wt_mat = repmat(wt_vec,max(wt_vec),1);
        [x_grid,y_grid] = meshgrid(1:size(wt_mat,2),1:size(wt_mat,1));
        x_ind_wt_mat(y_grid>wt_mat) = 0;
        y_ind_wt_mat(y_grid>wt_mat) = 0;
        x_ind_vec = reshape(x_ind_wt_mat,1,[]);
        y_ind_vec = reshape(y_ind_wt_mat,1,[]);
        x_ind_vec = x_ind_vec(x_ind_vec>0);
        y_ind_vec = y_ind_vec(y_ind_vec>0);                
        
        %%% conduct fitting
        y_vec_sp = (1:size(stripe_mat,1))';
        ap_center_guess = initial_centers(k,3); % initialize to cluster center
        [~, x_center_guess] = min(min(abs(fov_partitions(i).pixel_ap_id_mat-ap_center_guess))); % convert to pix
        x_diff_vec = abs(x_ind_vec-x_center_guess);
        x_ind_vec = x_ind_vec(x_diff_vec<x_tol);
        y_ind_vec = y_ind_vec(x_diff_vec<x_tol);
        pp = polyfit(y_ind_vec,x_ind_vec,3);
        x_poly = polyval(pp,y_vec_sp);        
        temp_center_mat(:,k,t_ind) = x_poly;
    end    
    % move backwards in time first
    orig_tp_vec = fliplr(1:t_critical/60-1); 
    for j = 1:length(orig_tp_vec)
        t_ind = find(t_fit_vec==orig_tp_vec(j));
        tp = orig_tp_vec(j)*60; 
        x_tol = 1.5*stripe_radius; 
        
        frame_array = zeros(size(stripe_mat)); 

        t_start = tp - t_window*60;
        t_stop = tp + t_window*60;

        t_fluo = fluo_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_x = xp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_y = yp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        idx = sub2ind(size(frame_array), t_y, t_x);
        frame_array(idx) = t_fluo;
        % apply smoothing kernel
        norm_ref_array = ones(size(frame_array));
        norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
        gauss_array = conv2(frame_array,g_kernel,'same');
        gauss_array = gauss_array./norm_ref_array;    

        %%% loop through stripes in set
        for k = 1:length(stripe_id_vec)
            stripe_id = stripe_id_vec(k);                          
            st_gauss =  gauss_array;
            st_gauss(stripe_mat~=stripe_id) = 0;   
            st_gauss_vec = reshape(st_gauss,1,[]);
            x_indices = repelem(1:size(stripe_mat,2),size(stripe_mat,1));                    
            y_indices = repmat(1:size(stripe_mat,1),1,size(stripe_mat,2));
            
            %%% apply weights
            wt_vec = floor(st_gauss_vec/100);
            if max(wt_vec) == 0
                warning('no fluorescence found. skipping')
                continue
            end
            x_ind_wt_mat = repmat(x_indices,max(wt_vec),1);
            y_ind_wt_mat = repmat(y_indices,max(wt_vec),1);
            wt_mat = repmat(wt_vec,max(wt_vec),1);
            [x_grid,y_grid] = meshgrid(1:size(wt_mat,2),1:size(wt_mat,1));
            x_ind_wt_mat(y_grid>wt_mat) = 0;
            y_ind_wt_mat(y_grid>wt_mat) = 0;
            x_ind_vec = reshape(x_ind_wt_mat,1,[]);
            y_ind_vec = reshape(y_ind_wt_mat,1,[]);
            x_ind_vec = x_ind_vec(x_ind_vec>0);
            y_ind_vec = y_ind_vec(y_ind_vec>0);                

            %%% conduct fitting
            y_vec_sp = (1:size(stripe_mat,1))';            
            x_fit_vec = [];
            y_fit_vec = [];
            for yp = 1:length(y_vec_sp)
                y = y_vec_sp(yp);
                xc = temp_center_mat(yp,k,t_ind+1);
                
                x_sub = x_ind_vec(y_ind_vec==y);
                y_sub = y_ind_vec(y_ind_vec==y);
                y_fit_vec = [y_fit_vec y_sub(abs(x_sub-xc)<=x_tol)];
                x_fit_vec = [x_fit_vec x_sub(abs(x_sub-xc)<=x_tol)];
            end            
            pp = polyfit(y_fit_vec,x_fit_vec,3);
            x_poly = polyval(pp,y_vec_sp);        
            temp_center_mat(:,k,t_ind) = x_poly;            
        end        
    end
    error('asfa')
    for i = 1
        if save_stripe_figs
            centroid_fig = figure;
            centroid_fig.Visible = 'off';
            centroid_fig.Position = [100 100 1024 256];
            hold on
            imagesc(gauss_array);
            scatter(reshape(center_mat,1,[]),repmat(1:size(stripe_mat,1),1,length(stripe_id_vec)),...
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
            saveas(centroid_fig,[flux_dynamics_path '/set_' num2str(i) '/ct_fluo_t' num2str(tp) '.tif'],'tif');
        end
        temp_fluo_array(:,:,j) = frame_array;
        temp_center_mat(:,:,j) = spline_mat;
        temp_disp_mat(:,:,j) = dispersion_mat;
    end    
    
    close all        
    %%%------------- Make time-dependent inference regions -------------%%%
    mean_fluo_mat = nanmean(temp_fluo_array,3);
    mean_center_mat = round(nanmean(temp_center_mat,3));
    stripe_id_mat_full = NaN(size(stripe_mat,1),size(stripe_mat,2),length(t_fit_vec));         
    for tp = 1:length(t_fit_vec)        
        center_mat = round(temp_center_mat(:,:,tp));
        for j = 1:length(stripe_id_vec)
            for k = 1:size(stripe_mat,1)
                stripe_id_mat_full(k,center_mat(k,j)-stripe_radius:...
                                center_mat(k,j)+stripe_radius,tp) = stripe_id_vec(j);
            end  
            if j == 1 %&& j ~= length(stripe_id_vec)
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,max(1,center_mat(k,j) - 3*stripe_radius):...
                                center_mat(k,j)-stripe_radius-1,tp) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                round(.5*(center_mat(k,j)+center_mat(k,j+1))),tp) =...
                                stripe_id_vec(j) + 1/3;                            
                end                                                 
            elseif j == length(stripe_id_vec)
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                center_mat(k,j)-stripe_radius-1,tp) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                min(1024,center_mat(k,j) + 3*stripe_radius),tp) =...
                                stripe_id_vec(j) + 1/3;                             
                end     
            else
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                    center_mat(k,j)-stripe_radius-1,tp) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                    round(.5*(center_mat(k,j)+center_mat(k,j+1))),tp) =...
                                    stripe_id_vec(j) + 1/3;                            
                end
            end
        end
    end
    mode_stripe_id_mat = mode(stripe_id_mat_full,3);
    mean_stripe_id_mat = mean(stripe_id_mat_full,3);
    stripe_pos_struct(i).stripe_id_mat = stripe_id_mat_full;
     
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
        xVec = xPos_vec_particle(particle_vec==ParticleID);
        yVec = yPos_vec_particle(particle_vec==ParticleID);
        t_trace = trace_struct(ind).time;
        f_trace = trace_struct(ind).fluo;
        t_trace = t_trace(~isnan(f_trace));
        tr_stripe_id_vec = zeros(1,length(t_trace))-1;
        for tp = 1:length(t_trace)
            xp = xVec(tp);
            yp = yVec(tp);
            filter = round(t_trace(tp)/60)==t_fit_vec;
            if sum(filter) == 0
                continue
            end
            s_id = stripe_id_mat_full(yp,xp,filter);
            if ~isnan(s_id)
                dispersion_factor = nanmean(temp_disp_mat(:,stripe_id_vec==round(s_id),filter));
                if dispersion_factor <= max_disp
                    tr_stripe_id_vec(tp) = s_id;
                else
                    tr_stripe_id_vec(tp) = -1;
                end
            else
                tr_stripe_id_vec(tp) = NaN;
            end
        end         
        fill_ind = find(-1==(tr_stripe_id_vec));
        id_vec = 1:length(tr_stripe_id_vec);
        bad_times = round(t_trace(fill_ind)/60);
        t_ref_vec = t_fit_vec;
        t_ref_vec(ismember(t_fit_vec,bad_times)) = Inf;
        if ~isempty(fill_ind) 
            % Assign early time points to nearest valid dynamic bin        
            for k = 1:length(fill_ind)
                [~, nn] = min(abs(t_ref_vec - round(t_trace(fill_ind(k))/60))); % find nearest valid point
                tr_stripe_id_vec(fill_ind(k)) = stripe_id_mat_full(yVec(fill_ind(k)),xVec(fill_ind(k)),nn);
            end            
        end        
        trace_struct(ind).stripe_id_inf = mode(tr_stripe_id_vec);
        trace_struct(ind).stripe_id_vec = tr_stripe_id_vec;
    end    
    
    %%% classify nuclei
    single_set_vec = [schnitz_struct.setID];
    set_indices = find(single_set_vec==i);
    
    for m = 1:length(set_indices)
        ind = set_indices(m);
        ncID = schnitz_struct(ind).ncID;
        xVec = round(xPos_vec_nc(nc_vec==ncID));
        yVec = round(yPos_vec_nc(nc_vec==ncID));
        t_trace = schnitz_struct(ind).time;                
        nc_stripe_id_vec = NaN(1,length(t_trace));
        for tp = 1:length(t_trace)
            xp = xVec(tp);
            yp = yVec(tp);
            filter = round(t_trace(tp)/60)==t_fit_vec;
            if sum(filter) == 0
                continue
            end
            s_id = stripe_id_mat_full(yp,xp,filter);
            if ~isnan(s_id)
                dispersion_factor = nanmean(temp_disp_mat(:,stripe_id_vec==round(s_id),filter));
                if dispersion_factor <= max_disp
                    nc_stripe_id_vec(tp) = s_id;
                else
                    nc_stripe_id_vec(tp) = -1;
                end
            else
                nc_stripe_id_vec(tp) = NaN;
            end
        end         
        fill_ind = find(-1==(nc_stripe_id_vec));
        id_vec = 1:length(nc_stripe_id_vec);
        bad_times = round(t_trace(fill_ind)/60);
        t_ref_vec = t_fit_vec;
        t_ref_vec(ismember(t_fit_vec,bad_times)) = Inf;
        if ~isempty(fill_ind) 
            % Assign early time points to nearest valid dynamic bin        
            for k = 1:length(fill_ind)
                [~, nn] = min(abs(t_ref_vec - round(t_trace(fill_ind(k))/60))); % find nearest valid point
                nc_stripe_id_vec(fill_ind(k)) = stripe_id_mat_full(yVec(fill_ind(k)),xVec(fill_ind(k)),nn);
            end            
        end                
        schnitz_struct(ind).stripe_id_inf = mode(nc_stripe_id_vec);
        schnitz_struct(ind).stripe_id_vec = nc_stripe_id_vec;
    end
    stripe_pos_struct(i).t_vec = t_fit_vec;
    disp(['Completed ' num2str(i) ' of ' num2str(length(set_index))])    
end
save(stripe_save_name, 'stripe_pos_struct')
save(trace_name,'trace_struct')
save(nucleus_name,'schnitz_struct')