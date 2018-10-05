% Script to make hi res parameter value maps
% Applies sliding window to soft-decode maps returned for each trace by HMM
% inference
%%%-----------------------Load Inference Results-------------------------%%

addpath('../utilities/');
clear 
close all
%%%%%%-----Set System Params
w = 7; %memory assumed for inference
K = 3; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
fluo_field = 1;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 30;
t_inf = 40;
% set write path 
DataPath = '../../dat/revisions/anti_corr/';
mkdir(DataPath)
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
% id variables
datatype = 'weka';
inference_type = 'dp';
project = 'eve7stripes_inf_2018_04_28'; %project identifier

%Generate filenames and writepath
id_string = [ '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins) ...
    '/K' num2str(K) '_t_window' num2str(round(t_window)) '_t_inf' num2str(round(t_inf)) '_' inference_type '/']; 

DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\eve7stripes_data\inference_out\';
% DropboxFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
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
load(['..\..\dat\' project '\inference_traces_' project '_dT20.mat']);

%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 1;
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'    
    load([folder_path filenames{f}]);
    if output.skip_flag == 1 
        continue
    end
    for fn = fieldnames(output)'
        glb_all(f_pass).(fn{1}) = output.(fn{1});
    end
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
% inf_stripe_vec = [];
mismatch_err = 0;
for i = 1:length(glb_all)
    traces = glb_all(i).traces;
    particles = glb_all(i).particle_ids;        
    inf_particle_id_vec = [inf_particle_id_vec particles];
%     inf_trace_id_vec = [inf_trace_id_vec NaN(1,length(particles))];
    inf_time_vec = [inf_time_vec repelem(glb_all(i).t_inf,length(particles))];
%     inf_stripe_vec = [inf_stripe_vec repelem(bin_range_vec(round(bin_map_vec,1)==glb_all(i).stripe_id),length(particles))];
    inf_id_vec = [inf_id_vec repelem(i,length(particles))];    
    sub_inf_id_vec = [sub_inf_id_vec 1:length(particles)];    
end

trace_particle_vec = [trace_struct_final.ParticleID];
inf_time_index = unique([glb_all.t_inf]);
inf_stripe_index = unique([glb_all.stripe_id]);
% params for soft decode inf
window_size = 3*60; % number of preceding and succeeding minutes to include in sliding window
hmm_window_size = glb_all(1).t_window/2;

details = struct;

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
%% Compile Trace-Specific Estimates for On and Off Rates
t_min = 25; % earliest time to use

transition_prob_mat = NaN(numel(trace_particle_vec),K^2);    
transition_rate_fit_mat = NaN(numel(trace_particle_vec),K^2);    
transition_rate_mat = NaN(numel(trace_particle_vec),K^2);    
tr_pt_id_vec = NaN(numel(trace_particle_vec),1);
tr_nc_id_vec = NaN(numel(trace_particle_vec),1);
tr_ap_vec = NaN(numel(trace_particle_vec),1);
tr_stripe_id_vec = NaN(numel(trace_particle_vec),1);
promoter_state_mat = NaN(numel(trace_particle_vec),K);    
fluo_vec = NaN(numel(trace_particle_vec),1);    
time_vec = NaN(numel(trace_particle_vec),1);    

parfor tr = 1:numel(trace_particle_vec)  
    ParticleID = trace_particle_vec(tr);
    % check to see if particle was included in inference runs
    pt_indices = find(inf_particle_id_vec==ParticleID);
    if isempty(pt_indices) % skip if not present
        continue
    end
    tt = trace_struct_final(trace_particle_vec==ParticleID).time_interp;
    t_filter = tt >= t_min*60;
    tr_nc_id_vec(tr) = trace_struct_final(tr).ncID;        
    tr_pt_id_vec(tr) = ParticleID;        
    tr_ap_vec(tr) = nanmean(trace_struct_final(tr).ap_vector_interp(t_filter));
    tr_stripe_id_vec(tr) = mode(trace_struct_final(tr).stripe_id_vec_interp(t_filter));
    fluo_vec(tr) = nanmean(trace_struct_final(tr).fluo_interp(t_filter));        
    % get particle indices    
    inf_ids = inf_id_vec(pt_indices);
    inf_sub_ids = sub_inf_id_vec(pt_indices);
    % take average of ss and s matrices for relevant time steps
    p_zz_tot = zeros(K,K);
    p_z_tot = zeros(K,1);                
    for id = 1:numel(inf_ids)
        [~, r_sort] = sort(glb_all(inf_ids(id)).r);
        p_zz = glb_all(inf_ids(id)).soft_struct.p_zz_log_soft{inf_sub_ids(id)};
        p_z = glb_all(inf_ids(id)).soft_struct.p_z_log_soft{inf_sub_ids(id)};        
        % time limits for inclusion              
        p_times = glb_all(inf_ids(id)).particle_times{inf_sub_ids(id)};
        pt_filter = p_times >=t_min*60;
        p_zz_tot = p_zz_tot + sum(exp(p_zz(r_sort,r_sort,pt_filter(2:end))),3);                    
        p_z_tot = p_z_tot + sum(exp(p_z(r_sort,pt_filter)),2);        
    end
    A_trace = p_zz_tot./sum(p_zz_tot);
    p_trace = p_z_tot / sum(p_z_tot(:));

    transition_prob_mat(tr,:) = reshape(A_trace,1,[]);
    promoter_state_mat(tr,:) = p_trace;
    % calculate transition and initiation rates
    R_trace = prob_to_rate(A_trace,Tres);                       
    %Check for imaginary and negative elements. If present, perform rate
    %fitting    
    R_out = R_trace;
    if ~isreal(R_out) || sum(R_out(:)<0)>K        
        warning('Performing Rate Fitting')
        out = prob_to_rate_fit_sym(A_trace, Tres, 'gen', .005, 1);            
        R_out = out.R_out;            
    end
    transition_rate_mat(tr,:) = reshape(R_trace,1,[]);
    transition_rate_fit_mat(tr,:) = reshape(R_out,1,[]);
    if mod(tr,50) == 0
        disp(['Completed ' num2str(tr) ' of ' num2str(numel(trace_struct_final))])
    end
end
%%
img_ft = imag(transition_rate_mat(:,2))==0 & imag(transition_rate_mat(:,4))==0;
nan_ft = ~isnan(transition_rate_mat(:,2));

kon_orig = real(transition_rate_mat(:,2)/2);
koff_orig = real(transition_rate_mat(:,4));

kon_fit = real(transition_rate_fit_mat(:,2)/2);
koff_fit = real(transition_rate_fit_mat(:,4));


out_table = array2table([tr_nc_id_vec(nan_ft) tr_ap_vec(nan_ft) img_ft(nan_ft) kon_orig(nan_ft) ...
    koff_orig(nan_ft) kon_fit(nan_ft) koff_fit(nan_ft) fluo_vec(nan_ft)],...
    'VariableNames', {'ncID','AP','PureReal','kon_orig','koff_orig','kon_fit','koff_fit','trace_fluo'});
writetable(out_table,[DataPath 'burst_rates_exp.csv'])

%%
close all
cm = jet(128);
check_fig = figure;
hold on
scatter(out_table.kon_orig,out_table.koff_orig,'MarkerFacecolor',cm(30,:),'MarkerFaceAlpha',.1,...
    'MarkerEdgeAlpha',0)
scatter(out_table.kon_orig(out_table.PureReal==1),out_table.koff_orig(out_table.PureReal==1),...
    'MarkerFacecolor',cm(120,:),'MarkerFaceAlpha',.1, 'MarkerEdgeAlpha',0)
legend('Fit Rates','Real Rates')
axis([0 .06 0 .06])