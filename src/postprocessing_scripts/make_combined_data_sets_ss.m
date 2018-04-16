% Script to generate combined data set containing inference traces alongside
% inference results
addpath('../utilities/');
clear 
close all
%%%%%%-----Set System Params
w = 7; %memory assumed for inference
K = 2; %states used for final inference
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
t_start_stripe = 25; 
%-----------------------------ID Variables--------------------------------%

% id variables
datatype = 'weka';
inference_type = 'dp';
project = 'eve7stripes_inf_2018_03_27_final'; %project identifier

%Generate filenames and writepath
id_thing = [ '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins)  '/']; 


DataPath = ['../../dat/' project '/' id_thing '/K' num2str(K) '_summary_stats/' ];
% load inference summary 
load([DataPath 'hmm_results_t_window' num2str(t_window) '_t_inf' num2str(t_inf)  '.mat'])
% load viterbi fits
load([DataPath 'viterbi_fits_t_window' num2str(t_window) '_t_inf' num2str(t_inf)  '.mat'])
v_specific = viterbi_fit_struct;
load([DataPath 'viterbi_fits_t_window50_t_inf25.mat'])
v_agg = viterbi_fit_struct;

% load traces
load(['..\..\dat\' project '\inference_traces_' project '_dT' num2str(Tres) '.mat']);
load(['..\..\dat\' project '\inference_nuclei_' project '_dT' num2str(Tres) '.mat']);
%%
%%% nucleus indexing vectors
nc_ncID_vec = [nuclei_clean.ncID];
nc_trID_vec = [nuclei_clean.ParticleID];
nc_set_vec = [nuclei_clean.setID];
%%% particle indexing vectors
tr_nc_vec = [trace_struct_final.ncID];
tr_particle_vec = [trace_struct_final.ParticleID];
tr_trID_vec = [trace_struct_final.ParticleID];
tr_set_vec = [trace_struct_final.setID];
%%% viterbi identifiers
va_particle_vec = [v_agg.ParticleID];
vs_particle_vec = [v_specific.ParticleID];

hmm_stripe_index = [hmm_results.binID]; % inf region ref vec
header = {'nucleus_id','particle_id', 'set_id', 'ap', 'xPos', 'yPos',...
          'stripe_id','inf_flag', 'time', 'fluo', 'tr_stripe_id', 'v_state_sp', 'v_fluo_sp'...
          'v_state_agg', 'v_fluo_agg', 'k_on', 'k_off', 'initiation_rate_off', 'initiation_rate_on'};

longform_data = NaN(length([nuclei_clean.time_interp]),length(header));
entry_index = 0;
for a = 1:length(nc_ncID_vec)  % use nuclei as basis for iteration    
    %%% extract nucleus details
    ncID = nc_ncID_vec(a);
    ParticleID = nc_trID_vec(a);
    setID = nc_set_vec(a);    
    nc_ap_vec = nuclei_clean(a).ap_vector_interp'; 
    nc_time_vec = nuclei_clean(a).time_interp'; 
    nc_xp_vec = nuclei_clean(a).xPos_interp'; 
    nc_yp_vec = nuclei_clean(a).yPos_interp';
    nc_stripe_vec = nuclei_clean(a).stripe_id_vec_interp';
    nc_setID_vec = repelem(setID,length(nc_stripe_vec))';
    nc_particle_vec = repelem(ParticleID,length(nc_stripe_vec))';
    nc_nucleus_vec = repelem(ncID,length(nc_stripe_vec))';
    nc_inf_flag_vec = repelem(0,length(nc_stripe_vec))';
    if ~isnan(ParticleID)
        % extract relevant particle metrics    
        trace_time = trace_struct_final(tr_trID_vec==ParticleID).time_interp;
        trace_fluo = trace_struct_final(tr_trID_vec==ParticleID).fluo_interp;                
        tr_inf_flag = trace_struct_final(tr_trID_vec==ParticleID).inference_flag;                
        tr_stripe_id_vec = trace_struct_final(tr_trID_vec==ParticleID).stripe_id_vec_interp(trace_time>t_start_stripe*60);    
        stripe_id = mode(tr_stripe_id_vec);  
        % longform id vectors        
        tr_stripe_id_long = repelem(stripe_id,length(trace_time))';    
        hmm_bin = hmm_results(round(hmm_stripe_index,1)==round(stripe_id,1));
        if ~isempty(hmm_bin)            
            R_mat = repmat(hmm_bin.R_fit_mean,length(trace_time),1);
            R_mat = R_mat(:,2:3)/60;
            r_mat = repmat(hmm_bin.initiation_mean'/60,length(trace_time),1);
            va_id = find(va_particle_vec==ParticleID);
            if isempty(va_id)
                va_mat = NaN(length(trace_time),2);
            else
                v_fit_a = v_agg(va_id).v_fit;
                va_mat = [double(v_fit_a.z_viterbi') v_fit_a.fluo_viterbi'];
            end
            vs_id = find(vs_particle_vec==ParticleID);
            if isempty(vs_id)
                vs_mat = NaN(length(trace_time),2);
            else
                v_fit_s = v_specific(vs_id).v_fit;
                vs_mat = [double(v_fit_s.z_viterbi') v_fit_s.fluo_viterbi'];
            end                
            % make longform ID vectors        
            trace_mat = [ trace_fluo'  double(tr_stripe_id_long)   double(vs_mat) double(va_mat) R_mat double(r_mat)];
        else
            trace_mat = [ trace_fluo'  tr_stripe_id_long   NaN(length(trace_fluo),8)];
        end
%         trace_mat = trace_mat(1:end-w-1,:);          
    end    
    nc_mat = [nc_nucleus_vec nc_particle_vec double(nc_setID_vec) nc_ap_vec double(nc_xp_vec)...
            double(nc_yp_vec) double(nc_stripe_vec) double(nc_inf_flag_vec) double(nc_time_vec)];
    if ~isnan(ParticleID)
        nc_tr_filter = ismember(round(nc_time_vec),round(trace_time));
%         nc_mat(nc_tr_filter,7) = tr_stripe_id_long;
        nc_mat(nc_tr_filter,8) = tr_inf_flag;
        nc_mat(~nc_tr_filter,2) = NaN;
    end    
    ind_vec = entry_index+1:entry_index+size(nc_mat,1); % align to main set    
    %%% add entries to master set
   longform_data(ind_vec,1:size(nc_mat,2)) = nc_mat;
   if ~isnan(ParticleID)
       longform_data(ind_vec(nc_tr_filter),size(nc_mat,2)+1:end) = trace_mat;
   end
   entry_index = ind_vec(end);
end
% remove NaN rows
% trace_mat_long = trace_mat_long(nan_index,:);
csvwrite_with_headers([DataPath '\eve_data_longform_w_nuclei.csv'], ...
                       longform_data, header,9); 
                 
save([DataPath 'eve_data_longform_w_nuclei.mat'],'longform_data')                   