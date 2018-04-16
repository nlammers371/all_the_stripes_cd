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
t_window = 15;
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = [bin_range_vec stripe_range(i) + j/3 - 2/3];
    end
end

% id variables
datatype = 'weka';
inference_type = 'set_bootstrap_results';
project = 'eve7stripes_inf_2018_02_20'; %project identifier

%Generate filenames and writepath
% truncated_inference_w7_t20_alpha14_f1_cl1_no_ends1
id_thing = [ '/truncated_inference_w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins) ...
    '/states' num2str(K) '/t_window' num2str(round(t_window)) '/' inference_type '/']; 

% DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\eve7stripes_data\inference_out\';
DropboxFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
folder_thing =  [DropboxFolder '/' project '/' id_thing];

OutPath = ['../../dat/' project '/' id_thing];
mkdir(OutPath)
load('..\..\dat\eve7stripes_inf_2018_02_20\w7_t20_alpha14_f1_cl1_no_ends1_tbins1\states2\t_window30\dp_bootstrap_results\viterbi_fits.mat')
load('..\..\dat\eve7stripes_inf_2018_02_20\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat');

viterbi_particle_vec = [viterbi_fit_struct.ParticleID];

%---------------------------------Read in Files---------------------------%
files = dir(folder_thing);
file_things = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
       ~isempty(strfind(files(i).name,['K' num2str(K)]))
        file_things = [file_things {files(i).name}];
    end
end

if isempty(file_things)
    error('No file with specified inference parameters found')
end
%%% load inference traces 


%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 0;

for f = 1:length(file_things)
    DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\eve7stripes_data\inference_out\';
    % DropboxFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
%     folder_thing =  [DropboxFolder '/' project '/' id_thing];
    % load the eve validation results into a structure array 'output'    
    load([folder_thing file_things{f}]);
    if output.skip_flag == 1 
        continue
    end
    for fn = fieldnames(output)'
        glb_all(f_pass).(fn{1}) = output.(fn{1});
    end
    glb_all(f_pass).source = file_things{f};        
    f_skip = f - f_pass
    f_pass = f_pass + 1    
end

% define indexing vectors
inf_id_vec = [];
sub_inf_id_vec = [];
inf_trace_id_vec = [];
inf_particle_id_vec = [];
inf_time_vec = [];

for i = 1:length(glb_all)
    traces = glb_all(i).traces;
    particles = glb_all(i).particle_ids;        
    inf_particle_id_vec = [inf_particle_id_vec particles];
    inf_time_vec = [inf_time_vec repelem(glb_all(i).t_inf,length(particles))];
    inf_id_vec = [inf_id_vec repelem(i,length(particles))];    
    sub_inf_id_vec = [sub_inf_id_vec 1:length(particles)];    
end
inf_particle_index = unique(inf_particle_id_vec);
trace_particle_vec = [trace_struct_final.ParticleID];
inf_time_index = unique([glb_all.t_inf]);
inf_stripe_index = unique([glb_all.stripe_id]);
% params for soft decode inf
hmm_window_size = glb_all(1).t_window/2;

for a = 1:length(inf_particle_index)            
    ParticleID = inf_particle_index(a);
    % extract relevant particle metrics
    setID = trace_struct_final(trace_particle_vec==ParticleID).setID;
    trace_time = trace_struct_final(trace_particle_vec==ParticleID).time_interp;
    trace_fluo = trace_struct_final(trace_particle_vec==ParticleID).fluo_interp;                
    stripe_id_vec = trace_struct_final(trace_particle_vec==ParticleID).stripe_id_vec_interp;        
    
    promoter_state_mat = NaN(length(trace_time),K);
    initiation_rate_mat = NaN(length(trace_time),K);
    transition_rate_mat = NaN(length(trace_time),K*(K-1));
    for i = 1:length(trace_time)     
        inf_time = trace_time(i);
        % inference results that fall into time window
        inf_filter = inf_time_vec-hmm_window_size<inf_time&...
            inf_time_vec+hmm_window_size>=inf_time;  
                
        p_filter = inf_filter&inf_particle_id_vec==ParticleID; % find inference runs with particle in current time window    
        inf_ids = inf_id_vec(p_filter);
        inf_sub_ids = sub_inf_id_vec(p_filter);
        % take average of ss and s matrices for relevant time steps        
        p_z_tot = zeros(K,1); 
        r_avg = zeros(K,1);
        R_avg = zeros(K*(K-1),1);
        for id = 1:length(inf_ids)
            [~, r_sort] = sort(glb_all(inf_ids(id)).r);
            A = glb_all(inf_ids(id)).A_mat;
            R = prob_to_rate(A,Tres);                        
            %Check for imaginary and negative elements. If present, perform rate
            %fit       
            R_out = R;
            if ~isreal(R) || sum(R(:)<0)>K                
                out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
                R_out = out.R_out;            
            end
            R_avg = R_avg + reshape(R_out(eye(K)~=1),[],1);
            r_avg = r_avg + glb_all(inf_ids(id)).r;
            p_z = glb_all(inf_ids(id)).soft_struct.p_z_log_soft{inf_sub_ids(id)};
            p_times = glb_all(inf_ids(id)).particle_times{inf_sub_ids(id)};            
            pt_filter = p_times == trace_time(i);            
            p_z_tot = p_z_tot + sum(exp(p_z(r_sort,pt_filter)),2);            
        end        
        p_slice = p_z_tot / sum(p_z_tot(:));
        r_slice = r_avg / sum(p_z_tot(:));
        R_slice = R_avg / sum(p_z_tot(:));              
        promoter_state_mat(i,:) = p_slice;
        transition_rate_mat(i,:) = R_slice; % record rate
        initiation_rate_mat(i,:) = r_slice;                
    end
    pID_vec = repelem(ParticleID,length(trace_time))';
    setID_vec = repelem(setID,length(trace_time))';
    trace_mat = [pID_vec setID_vec stripe_id_vec' trace_time' trace_fluo' ...
        promoter_state_mat  transition_rate_mat initiation_rate_mat];
    trace_mat = trace_mat(1:end-w-1,:);
    if a == 1
        trace_mat_long = trace_mat;
    else
        trace_mat_long = [trace_mat_long ; trace_mat];
    end
end
%%
% remove NaN rows
nan_index = ~isnan(sum(trace_mat_long,2));
trace_mat_long = trace_mat_long(nan_index,:);
header = {'particle_id', 'set_id', 'stripe_id','time', 'fluo', 'p_off', 'p_on'...
          , 'k_on', 'k_off', 'initiation_rate_off', 'initiation_rate_on'};
csvwrite_with_headers([OutPath '\eve_data_longform.csv'], ...
                       trace_mat_long, header); 
save([OutPath 'eve_data_longform.mat'],'trace_mat_long')                   