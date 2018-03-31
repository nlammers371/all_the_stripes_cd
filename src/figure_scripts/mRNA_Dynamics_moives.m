% Script to Make Movies Depicting Instantaneous Stripe Dynamics and mRNA
% Accumulation

close all
clear 

% set filenames
project = 'eve7stripes_inf_2018_03_27'; %Project Identifier

fig_path = ['../../fig/experimental_system/' project '/stripe_dynamics/'];
data_path = ['../../dat/' project '/']; % data mat directory

flux_dynamics_path = [fig_path 'stripe_mRNA_dynamics/'];
mkdir(flux_dynamics_path);
trace_name = [data_path 'inference_traces_' project '_dT20.mat']; % names for compiled trace struct
nucleus_name = [data_path 'inference_nuclei_' project '_dT20.mat']; % names for compiled elipse struct
fov_name = [data_path 'fov_partitions_' project '.mat'];
stripe_name = [data_path 'stripe_pos_' project '.mat'];

% load processed datasets
load(trace_name); % particle info
load(fov_name); % ap and stripe info at pixel level
load(nucleus_name);
load(stripe_name)

% FOV size
xDim = 1024;
yDim = 256;
[x_ref_mat,y_ref_mat] = meshgrid(1:xDim,1:yDim);
%% stripe dynamics movies
cm = jet(128);
increment = floor(size(cm,1)/7);
%Array to store color mappings
stripe_colors = zeros(7, 3);
for i = 1:7
    stripe_colors(i,:) = cm(1+(i-1)*increment,:);
end
%%% time params
frame_rate = 1; %minutes
half_life = 7; % mRNA half life (minutes)
t_vec = 5:frame_rate:50; % minutes
decay_prob = frame_rate/half_life; % likelihood of decay from frame to frame
%%% contrast params
MaxFluo = 3500;
MaxmRNA = 1e7;
% nucleus_radius = 25; % pixels 
MaxRadius = 35;
%%% convolution kernel to spread spot
g_kernel_radius = 4; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 2; % this is kind of arbitrary
[x_ref, y_ref] = meshgrid(1:2*g_kernel_radius+1,1:2*g_kernel_radius+1);
x_ref = x_ref - g_kernel_radius - 1;
y_ref = y_ref - g_kernel_radius - 1;
r_mat = sqrt(x_ref.^2 + y_ref.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel
g_kernel(r_mat>g_kernel_radius) = 0;
% make indexing vectors for traces
xPos_vec_particle = [];
yPos_vec_particle = [];
set_vec_particle = [];
fluo_vec_particle = [];
time_vec_particle = [];
stripe_vec_particle = [];
nc_vec_particle = [];
particle_vec = [];
for i = 1:length(trace_struct_final)    
    particle_vec = [particle_vec repelem(trace_struct_final(i).ParticleID,length(trace_struct_final(i).yPos_interp))];
    nc_vec_particle = [nc_vec_particle repelem(trace_struct_final(i).ncID,length(trace_struct_final(i).yPos_interp))];
    fluo_vec_particle = [fluo_vec_particle trace_struct_final(i).fluo_interp];
    time_vec_particle = [time_vec_particle trace_struct_final(i).time_interp];
    xPos_vec_particle = [xPos_vec_particle trace_struct_final(i).xPos_interp];
    yPos_vec_particle = [yPos_vec_particle trace_struct_final(i).yPos_interp];    
    stripe_vec_particle = [stripe_vec_particle  trace_struct_final(i).stripe_id_vec_interp];    
    set_vec_particle = [set_vec_particle repelem(trace_struct_final(i).setID, length(trace_struct_final(i).yPos))];    
end
% same for nuclei
xPos_vec_nc = [];
yPos_vec_nc = [];
stripe_vec_nc = [];
set_vec_nc = [];
fluo_vec_nc = [];
time_vec_nc = [];
nc_vec = [];
for i = 1:length(nuclei_clean)        
    nc_vec = [nc_vec repelem(nuclei_clean(i).ncID,length(nuclei_clean(i).yPos_interp))];
    time_vec_nc = [time_vec_nc nuclei_clean(i).time_interp];
    stripe_vec_nc = [stripe_vec_nc nuclei_clean(i).stripe_id_vec_interp];
    xPos_vec_nc = [xPos_vec_nc nuclei_clean(i).xPos_interp];
    yPos_vec_nc = [yPos_vec_nc nuclei_clean(i).yPos_interp];    
    set_vec_nc = [set_vec_nc repelem(nuclei_clean(i).setID, length(nuclei_clean(i).yPos_interp))];    
end
set_index = unique(set_vec_particle);
%%
stripe_dynamics_struct = struct; % save stripe location arrays
for i = 2:length(set_index)
    stripe_id_array = stripe_pos_struct(i).stripe_id_mat;
    % particle info
    tr_xp_set_vec = xPos_vec_particle(set_vec_particle==set_index(i));
    tr_yp_set_vec = yPos_vec_particle(set_vec_particle==set_index(i));
    tr_fluo_set_vec = fluo_vec_particle(set_vec_particle==set_index(i));
    tr_time_set_vec = time_vec_particle(set_vec_particle==set_index(i));
    tr_pid_set_vec = particle_vec(set_vec_particle==set_index(i));
    tr_nid_set_vec = nc_vec_particle(set_vec_particle==set_index(i));
    % nc info
    nc_xp_set_vec = xPos_vec_nc(set_vec_nc==set_index(i));
    nc_yp_set_vec = yPos_vec_nc(set_vec_nc==set_index(i));    
    nc_time_set_vec = time_vec_nc(set_vec_nc==set_index(i));
    nc_id_set_vec = nc_vec(set_vec_nc==set_index(i));
    ap_mat = fov_partitions(i).pixel_ap_id_mat;    
    %%% handle inversions
    mean_ap_vec = nanmean(ap_mat);
    if mean_ap_vec(1) > mean_ap_vec(end)                
        tr_xp_set_vec = xDim - tr_xp_set_vec + 1;
        nc_xp_set_vec = xDim - nc_xp_set_vec + 1;        
    end    
    set_t_vec = stripe_pos_struct(i).t_vec; % time vector corresponding to stripe array
    temp_fluo_array = zeros(yDim,xDim,length(t_vec)); % store instantaneous Fluo
    temp_mRNA_array = zeros(yDim,xDim,length(t_vec)); %store accumulated mRNA    
    for CurrentFrame = 1:length(t_vec)
        ind = max(1,t_vec(CurrentFrame)-set_t_vec(1) + 1);
        stripe_id_mat = stripe_id_array(:,:,ind);
        frame_array = zeros(yDim,xDim);
        t = t_vec(CurrentFrame)*60;
        t_start = t - 60;
        t_stop = t;
        % get set of values for current time period
        t_filter_tr = tr_time_set_vec>=t_start&tr_time_set_vec<t_stop;
        t_filter_nc = nc_time_set_vec>=t_start&nc_time_set_vec<t_stop;
        tr_t_fluo = tr_fluo_set_vec(t_filter_tr);
        tr_t_x = tr_xp_set_vec(t_filter_tr);
        tr_t_y = tr_yp_set_vec(t_filter_tr);
        tr_t_pid = tr_pid_set_vec(t_filter_tr);
        tr_t_nid = tr_nid_set_vec(t_filter_tr);
        
        nc_t_id = nc_id_set_vec(t_filter_nc);
        nc_t_y = nc_yp_set_vec(t_filter_nc);
        nc_t_x = nc_xp_set_vec(t_filter_nc);
        
        nc_index = unique(nc_t_id); % unique list of extant nuclei
        particle_index = unique(tr_t_pid);
        for tr = 1:length(nc_index)       
            tr_filter = tr_t_nid==nc_index(tr);
            nc_filter = nc_t_id==nc_index(tr);
%             if length(unique(tr_t_pid(tr_filter))) > 1
%                 error('nc bs')
%             end
            xp = round(mean(nc_t_x(nc_filter)));
            yp = round(mean(nc_t_y(nc_filter)));
            fp = sum(tr_t_fluo(tr_filter));
            n_particle = length(unique(tr_t_pid(tr_t_nid==nc_index(tr))));            
            if ~isnan(xp)
                frame_array(yp,xp) = fp;
            end
        end
%         tr_idx = sub2ind(size(frame_array), tr_t_y, tr_t_x);
%         frame_array(round(tr_idx)) = tr_t_fluo;        
        % apply smoothing kernel
        norm_ref_array = ones(size(frame_array));
        norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
        gauss_array = conv2(frame_array,g_kernel,'same');
        gauss_array = gauss_array./norm_ref_array;
        temp_fluo_array(:,:,CurrentFrame) = gauss_array; % store instantaneous fluo
        % generate nuclear patches
        % find nucleu centers
        nc_x = zeros(1,length(nc_index));
        nc_y = zeros(1,length(nc_index));
        for nc = 1:length(nc_index)
            nc_x(nc) = round(mean(nc_t_x(nc_t_id==nc_index(nc))));
            nc_y(nc) = round(mean(nc_t_y(nc_t_id==nc_index(nc))));
        end        
        NucleusmRNAMat = zeros(yDim,xDim,3);        
        NucleusIDMat = zeros(yDim,xDim);
        NucleusDistMat = ones(yDim,xDim)*MaxRadius;    
        % Loop through ALL extant nuclei        
        for s = 1:length(nc_index)
            x =  nc_x(s);
            y =  nc_y(s);        
            distances = ((x_ref_mat-x).^2 + (y_ref_mat-y).^2).^.5;
            candidate_indices = NucleusDistMat > distances; 
            %assign IDs    
            NucleusIDMat(candidate_indices) = nc_index(s);
            NucleusDistMat(candidate_indices) = distances(candidate_indices);            
        end
        %%% assign mRNA values to each nucleus
        new_mRNA = zeros(size(NucleusIDMat));
        f_slice = zeros(size(NucleusIDMat));
        if CurrentFrame > 1
            f_slice = temp_fluo_array(:,:,CurrentFrame-1);
            for nc = 1:length(nc_index)
                nc_id = nc_index(nc);
                nc_idx = find(NucleusIDMat==nc_id);                
                new_mRNA(nc_idx) = sum(f_slice(nc_idx)); % whole nuclear region gets "credit" for production
            end
            mRNA_slice = temp_mRNA_array(:,:,CurrentFrame-1);
            temp_mRNA_array(:,:,CurrentFrame) = (1-decay_prob)*mRNA_slice + new_mRNA; 
        end
        %%% assign colors to nucleus patches
        for nc = 1:length(nc_index)
            nc_id = nc_index(nc);
            mRNA_slice = temp_mRNA_array(:,:,CurrentFrame);
            nc_mRNA = mean(mRNA_slice(NucleusIDMat==nc_id));
            stripe_id = stripe_id_mat(nc_y(nc),nc_x(nc));
            frac = min(1,nc_mRNA/MaxmRNA);
            if isnan(stripe_id)
                stripe_color = .4*(1-frac) + [0 0 0] + [.95 .95 .95]*frac;
            else
                stripe_color = .4*(1-frac)+[0 0 0] + stripe_colors(round(stripe_id),:)*frac;
            end            
            %Record Fluorescence
            for k = 1:3
                filter = NucleusIDMat==nc_id;        
                nc_slice = NucleusmRNAMat(:,:,k);
                nc_slice(filter) = stripe_color(k);
                NucleusmRNAMat(:,:,k) = nc_slice;        
            end                
        end        
        %Now Draw Nucleus Borders for active nuclei        
%         window = 1;
        for nc = 1:length(nc_index)
            nc_id = nc_index(nc);
            stripe_id = stripe_id_mat(nc_y(nc),nc_x(nc));
            if isnan(stripe_id) || round(stripe_id)~=stripe_id
                stripe_color_edge = [.37 .37 .37];
                stripe_color_pt = [0 0 0] + .95;
                if ~isnan(stripe_id)
                    stripe_color_pt = stripe_colors(round(stripe_id),:);
                end
            else
                stripe_color_edge = stripe_colors(round(stripe_id),:);
                stripe_color_pt = stripe_colors(round(stripe_id),:);
            end                                     
            
            %get coordinates of nucleus patch
            b_kernel = ones(3,3);
            binary_map = NucleusIDMat==nc_id;
            norm_ref_array = ones(size(binary_map));
            norm_ref_array = conv2(norm_ref_array,b_kernel,'same');
            edge_array = conv2(binary_map,b_kernel,'same');
            edge_array = edge_array./norm_ref_array;
            edge_filter = edge_array < 1 & edge_array > .5; 
            spot_filter = NucleusIDMat==nc_id;
            spot_mask = gauss_array(NucleusIDMat==nc_id)/MaxFluo;
            for k = 1:3
                nc_slice = NucleusmRNAMat(:,:,k); 
                nc_orig = nc_slice(spot_filter);
                nc_slice(spot_filter) = stripe_color_pt(k);
                nc_slice(spot_filter) = nc_slice(spot_filter).*spot_mask + (1-spot_mask).*nc_orig;
                if ~isnan(max(stripe_color_edge))
                    nc_slice(edge_filter) = stripe_color_edge(k);
                end
                NucleusmRNAMat(:,:,k) = nc_slice;
%                 if sum(spot_filter(:)) ~= 0 && nc > 11 && k == 3
%                     error('asfsa')
%                 end
            end            
        end
        mRNA_fig = figure;
        imshow(NucleusmRNAMat) 
        title(['mRNA Dynamics Set ' num2str(i) ' (' num2str(t_vec(CurrentFrame)) ' min)'])
%         axis([0 yDim 0 xDim])
%         text(10,25,[iIndex(round(t/60),2),...
%         ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])
        ax = gca;
        ax.Visible = 'off';
        mkdir([flux_dynamics_path '/set_' num2str(i) '/'])
        saveas(mRNA_fig,[flux_dynamics_path '/set_' num2str(i) '/mRNA_mask_min' ...
            num2str(t_vec(CurrentFrame)) '.tif'],'tif');        
    end
end    
close all