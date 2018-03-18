% Script to make hi res parameter value maps
% Applies sliding window to soft-decode maps returned for each trace by HMM
% inference
%%%-----------------------Load Inference Results-------------------------%%

addpath('../utilities/');
clear all
close all
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
% DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\hmmm_data\inference_out\';
DropboxFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
folder_path =  [DropboxFolder '/' project '/' id_string];
OutPath = ['../../dat/' project '/' id_string];
FigPath = ['../../fig/experimental_system/' project '/' id_string];
mkdir(OutPath)
mkdir(FigPath)
%---------------------------------Read in Files---------------------------%
files = dir(folder_path);
filenames = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
       ~isempty(strfind(files(i).name,['K' num2str(K)]))
        filenames = [filenames {files(i).name}];
    end
end

if isempty(filenames)
    error('No file with specified inference parameters found')
end
%%% load inference traces 
load('E:\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_02_20\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat');

%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 1;
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'    
    load([folder_path filenames{f}]);
    if output.skip_flag == 1 
        continue
    elseif length(output.stripe_id) > 1
        continue
    elseif output.t_window ~= 900
        continue
    end
    for fn = fieldnames(output)'
        glb_all(f_pass).(fn{1}) = output.(fn{1});
    end
%     if glb_all(f_pass).stripe_id == round(glb_all(f_pass).stripe_id) - .1
%         glb_all(f_pass).stripe_id = round(glb_all(f_pass).stripe_id) - 1/3;
%     elseif glb_all(f_pass).stripe_id == round(glb_all(f_pass).stripe_id) + .1
%         glb_all(f_pass).stripe_id = round(glb_all(f_pass).stripe_id) + 1/3;
%     end
    glb_all(f_pass).source = filenames{f};        
    f_pass = f_pass + 1
%     catch
%         warning(['File ' num2str(f) ' failed to load'])
%     end
end
%%
% make a list of indices and particle ids contained in each inf entry
inf_id_vec = [];
sub_inf_id_vec = [];
inf_trace_id_vec = [];
inf_particle_id_vec = [];
inf_time_vec = [];
inf_stripe_vec = [];
mismatch_err = 0;
for i = 1:length(glb_all)
    traces = glb_all(i).traces;
    particles = glb_all(i).particle_ids;        
    inf_particle_id_vec = [inf_particle_id_vec particles];
%     inf_trace_id_vec = [inf_trace_id_vec NaN(1,length(particles))];
    inf_time_vec = [inf_time_vec repelem(glb_all(i).t_inf,length(particles))];
    inf_stripe_vec = [inf_stripe_vec repelem(bin_range_vec(round(bin_map_vec,1)==glb_all(i).stripe_id),length(particles))];
    inf_id_vec = [inf_id_vec repelem(i,length(particles))];    
    sub_inf_id_vec = [sub_inf_id_vec 1:length(particles)];    
end

trace_particle_vec = [trace_struct_final.ParticleID];
inf_time_index = unique([glb_all.t_inf]);
inf_stripe_index = unique([glb_all.stripe_id]);
% params for soft decode inf
tDelta = 60; % temporal step size (in minutes) 
sd_stripe_regions = (6*min(bin_range_vec):1:6*max(bin_range_vec))/6;
sd_time_regions = min(inf_time_index):tDelta:max(inf_time_index);
window_size = 3*60; % number of preceding and succeeding minutes to include in sliding window
hmm_window_size = glb_all(1).t_window/2;

details = struct;
min_dp_per_point_vec = repelem(150,length(sd_stripe_regions));
n_boots = 1;
% Make alpha correction kernel
alpha_kernel = [];
times = 0:Tres:w*Tres;
for n = 1:w
    t1 = times(n);
    t2 = times(n+1);
    alpha_kernel = [alpha_kernel ms2_loading_coeff_integral(alpha, w, Tres, t1, t2)];
end
alpha_kernel = fliplr(alpha_kernel)/Tres;
glb_index_vec = 1:length(glb_all);

transition_prob_mat= NaN(length(sd_time_regions),length(sd_stripe_regions),K^2);    
transition_rate_mat= NaN(length(sd_time_regions),length(sd_stripe_regions),K^2);    
promoter_state_mat = NaN(length(sd_time_regions),length(sd_stripe_regions),K);    
initiation_rate_mat = NaN(length(sd_time_regions),length(sd_stripe_regions),K);    
obs_mat_all = NaN(length(sd_time_regions),length(sd_stripe_regions));   
fluo_mat_all = NaN(length(sd_time_regions),length(sd_stripe_regions));    

for a = 1:length(sd_stripe_regions)            
    analysis_stripe_region = sd_stripe_regions(a);                 
    for b = 1:n_boots        
        % make arrays to store results            
        fluo_profile = zeros(1,length(sd_time_regions));
        fluo_counts = zeros(1,length(sd_time_regions));                              
        % loop through time points. select traces on the fly, allowing for
        % fact that trace region assignment can vary in time
        for i = 1:length(sd_time_regions)     
            inf_time = sd_time_regions(i);
            % find elligible inference results
            filter = inf_time_vec-hmm_window_size<inf_time&...
                inf_time_vec+hmm_window_size>inf_time &  inf_stripe_vec<=analysis_stripe_region+1/6 &...
                inf_stripe_vec>=analysis_stripe_region-1/6;            
            eligible_particles = unique(inf_particle_id_vec(filter));
            inf_particle_vec = [];            
            for m = 1:length(eligible_particles) % determine which particles are in appt sub region
                tt = trace_struct_final(trace_particle_vec==eligible_particles(m)).time;
                ft = trace_struct_final(trace_particle_vec==eligible_particles(m)).fluo;                
                stripe_id_vec = trace_struct_final(trace_particle_vec==eligible_particles(m)).stripe_id_vec;
                tt = tt(~isnan(ft));                
                mean_stripe_id = mean(stripe_id_vec(tt>=inf_time-window_size*60&tt<inf_time+window_size*60));
                if mean_stripe_id<analysis_stripe_region+1/6&&mean_stripe_id>=analysis_stripe_region-1/6
                    inf_particle_vec = [inf_particle_vec eligible_particles(m)];
                end
            end
            boot_vec = 1:length(inf_particle_vec);
            boot_ids = inf_particle_vec;%randsample(boot_vec,length(inf_particle_vec),true);
            % loop through filtered particle list. Randomly choose a
            % bootstrap from each relvant time window. Take average across
            % randomly selected bootstraps
            slice_zz = zeros(K,K);
            slice_z = zeros(K,1);
            slice_fluo_ct = zeros(length(inf_particle_vec),K);
            slice_fluo = zeros(length(inf_particle_vec),1);
            for p = 1:length(boot_vec)
                ParticleID = boot_ids(p);
                p_filter = filter&inf_particle_id_vec==ParticleID;          
                inf_ids = inf_id_vec(p_filter);
                inf_sub_ids = sub_inf_id_vec(p_filter);
                % take average of ss and s matrices for relevant time steps
                p_zz_tot = zeros(K,K);
                p_z_tot = zeros(K,1);
                fluo_tot = 0;                
                for id = 1:length(inf_ids)
                    [~, r_sort] = sort(glb_all(inf_ids(id)).r);
                    p_zz = glb_all(inf_ids(id)).soft_struct.p_zz_log_soft{inf_sub_ids(id)};
                    p_z = glb_all(inf_ids(id)).soft_struct.p_z_log_soft{inf_sub_ids(id)};
                    p_times = glb_all(inf_ids(id)).particle_times{inf_sub_ids(id)};
                    p_fluo = glb_all(inf_ids(id)).traces{inf_sub_ids(id)};
                    % time limits for inclusion
                    sd_start = inf_time-window_size;
                    sd_stop = inf_time+window_size;                     
%                     sd_filter = sum(); % usable trace times                    
                    pt_filter = p_times >=sd_start & p_times < sd_stop;
                    p_zz_tot = p_zz_tot + sum(exp(p_zz(r_sort,r_sort,pt_filter(2:end))),3);                    
                    p_z_tot = p_z_tot + sum(exp(p_z(r_sort,pt_filter)),2);
                    fluo_tot = fluo_tot + sum(p_fluo(pt_filter));                    
                end
                A_slice = p_zz_tot./sum(p_zz_tot);
                if ~isnan(max(A_slice(:)))
                    slice_zz = slice_zz + A_slice;
                end
                p_slice = p_z_tot / sum(p_z_tot(:));
                if ~isnan(max(p_slice(:)))
                    slice_z = slice_z + p_z_tot / sum(p_z_tot(:));
                    slice_fluo_ct(p,:) = p_z_tot / sum(p_z_tot(:));
                    slice_fluo(p) = fluo_tot / sum(p_z_tot(:));
                end                              
            end
            A = slice_zz ./ sum(slice_zz);
            if isnan(max(A(:)))
                continue
            end
            transition_prob_mat(i,a,:) = reshape(A,1,[]);
            promoter_state_mat(i,a,:) = slice_z ./ sum(slice_z);
            fluo_mat_all(i,a) = mean(slice_fluo);
            obs_mat_all(i,a) = sum(slice_z);
            % calculate transition and initiation rates
            R = prob_to_rate(A,Tres);                        
            %Check for imaginary and negative elements. If present, perform rate
            %fit       
            R_out = R;
            if ~isreal(R) || sum(R(:)<0)>K        
%                 disp('Pathological rate matrix found. Performing fitting')
                out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
                R_out = out.R_out;            
            end
            transition_rate_mat(i,a,:) = reshape(R_out,1,[]); % record rate
            % Emission values
            F_square = zeros(K,K);
            b_agg = zeros(1,K);
            if ~isempty(slice_fluo_ct)
                for k = 1:K
                    F_square(k,:) = sum(slice_fluo_ct.*repmat(slice_fluo_ct(:,k),1,K));
                    b_agg(k) = sum(slice_fluo.*slice_fluo_ct(:,k));
                end
                % solve system                
                initiation_rate_mat(i,a,:) = linsolve(F_square,b_agg')'; 
            end
        end        
    end
    disp(['Completed ' num2str(a) ' of ' num2str(length(sd_stripe_regions))])
    details(a).min_dp = min_dp_per_point_vec(a);       
    details(a).stripe_region = analysis_stripe_region;
    details(a).time_vec = sd_time_regions;
end
soft_decode_params.details = details;
soft_decode_params.initiation_rate_mat = initiation_rate_mat;
soft_decode_params.promoter_state_mat = promoter_state_mat;
soft_decode_params.transition_rate_mat = transition_rate_mat;
soft_decode_params.transition_prob_mat = transition_prob_mat;
soft_decode_params.stripe_region_vec = sd_stripe_regions;
soft_decode_params.mean_fluo_mat = fluo_mat_all;
save([OutPath '/soft_decode_params.mat'],'soft_decode_params')