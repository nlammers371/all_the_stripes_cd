% Script to generate combined data set containing inference traces alongside
% inference results
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
project = 'eve7stripes_inf_2018_04_28'; %project identifier

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
% load([DataPath 'viterbi_fits_t_window' num2str(t_window) '_t_inf' num2str(t_inf)  '_agg.mat'])
v_agg = viterbi_fit_struct;

% load traces
load(['..\..\dat\' project '\inference_traces_' project '_dT' num2str(Tres) '.mat']);
load(['..\..\dat\' project '\inference_nuclei_' project '_dT' num2str(Tres) '.mat']);

%%% nucleus indexing vectors
nc_ncID_vec = [nucleus_struct_final.ncID];
nc_trID_vec = [nucleus_struct_final.ParticleID];
nc_set_vec = [nucleus_struct_final.setID];
%%% particle indexing vectors
tr_nc_vec = [trace_struct_final.ncID];
tr_particle_vec = [trace_struct_final.ParticleID];
tr_trID_vec = [trace_struct_final.ParticleID];
tr_set_vec = [trace_struct_final.setID];
%%% viterbi identifiers
va_particle_vec = [v_agg.ParticleID];
vs_particle_vec = [v_specific.ParticleID];

hmm_stripe_index = [hmm_results.binID]; % inf region ref vec
header = {'nucleus_id','particle_id', 'set_id', 'ap_raw', 'ap_registered', 'xPos', 'yPos',...
          'stripe_id','inf_flag', 'time', 'fluo', 'tr_stripe_id', 'v_state_sp', 'v_fluo_sp'...
          'v_state_agg', 'v_fluo_agg', 'k_on', 'k_off', 'initiation_rate_off', 'initiation_rate_on'};

longform_data = NaN(length([nucleus_struct_final.time_interp]),length(header));
entry_index = 0;
for a = 1:length(nc_ncID_vec)  % use nuclei as basis for iteration    
    %%% extract nucleus details
    ncID = nc_ncID_vec(a);
    ParticleID = nc_trID_vec(a);
    setID = nc_set_vec(a);    
    nc_ap_vec_raw = nucleus_struct_final(a).ap_vector_interp'; 
    nc_ap_vec_registered = nucleus_struct_final(a).ap_vector_warped'; 
    nc_time_vec = nucleus_struct_final(a).time_interp'; 
    nc_xp_vec = nucleus_struct_final(a).xPos_interp'; 
    nc_yp_vec = nucleus_struct_final(a).yPos_interp';
    nc_stripe_vec = nucleus_struct_final(a).stripe_id_vec_interp';
    nc_setID_vec = repelem(setID,length(nc_stripe_vec))';
    nc_particle_vec = repelem(ParticleID,length(nc_stripe_vec))';
    nc_nucleus_vec = repelem(ncID,length(nc_stripe_vec))';
    nc_inf_flag_vec = repelem(NaN,length(nc_stripe_vec))';
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
            R_mat_full = repmat(hmm_bin.R_fit_mean,length(trace_time),1);
            if K == 2
                R_mat = R_mat_full(:,2:3)/60;
                r_mat = repmat(hmm_bin.initiation_mean'/60,length(trace_time),1);
            elseif K == 3
                effective_off = ((R_mat_full(:,4)+R_mat_full(:,6))./R_mat_full(:,4).*-R_mat_full(:,5).^-1 + ...
                                R_mat_full(:,6) ./ R_mat_full(:,4).*-R_mat_full(:,9).^-1).^-1;
                R_mat = [R_mat_full(:,2)/60 effective_off/60];                
                r_mat_full = repmat(hmm_bin.initiation_mean'/60,length(trace_time),1);
                s2 = -R_mat_full(1,9)^-1;
                s1 = -R_mat_full(1,5)^-1;
                eff_init = s1/(s1+s2)*r_mat_full(:,2)+s2/(s1+s2)*r_mat_full(:,3);
                r_mat = [r_mat_full(:,1) eff_init];                
            end
            
            va_id = find(va_particle_vec==ParticleID);
            if isempty(v_agg(va_id).v_fit)
                va_mat = NaN(length(trace_time),2);
            else
                v_fit_a = v_agg(va_id).v_fit;
                va_mat = [double(v_fit_a.z_viterbi') v_fit_a.fluo_viterbi'];
            end
            vs_id = find(vs_particle_vec==ParticleID);
            if isempty(v_specific(vs_id).v_fit)
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
    nc_mat = [nc_nucleus_vec nc_particle_vec double(nc_setID_vec) nc_ap_vec_raw nc_ap_vec_registered double(nc_xp_vec)...
            double(nc_yp_vec) double(nc_stripe_vec) double(nc_inf_flag_vec) double(nc_time_vec)];
    if ~isnan(ParticleID)
        nc_tr_filter = ismember(round(nc_time_vec),round(trace_time));
%         nc_mat(nc_tr_filter,7) = tr_stripe_id_long;
        nc_mat(nc_tr_filter,9) = tr_inf_flag;
%         nc_mat(~nc_tr_filter,9) = NaN;
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
%% Perform QC Checks
cm = jet(128);
% make sure inference and nucleus stripe id assignements are consistent
stripe_id_fig = figure;
scatter(longform_data(:,8), longform_data(:,12),40,'MarkerFaceColor',cm(40,:),'MarkerEdgeColor','black');
xlabel('nucleus stripe ID')
ylabel('trace inference stripe ID')
title('Consisitency of Stripe ID Variables')
saveas(stripe_id_fig,[DataPath '\stripe_consitency_check.png'],'png')

% check to see if dynamic bins track with activity over time
plot_times = 1:50;
nc_stripe_id_vec = longform_data(:,8);
tr_fluo_vec = longform_data(:,11);
tr_fluo_vec(isnan(tr_fluo_vec)) = 0;
v_state_vec = longform_data(:,13);

time_vec = round(longform_data(:,10)/60);
stripe_index = unique(nc_stripe_id_vec);
stripe_index = stripe_index(~isnan(stripe_index));
%%% track activity over time
f_activity_mat = NaN(length(plot_times),length(stripe_index));
v_activity_mat = NaN(length(plot_times),length(stripe_index));
for t = 1:length(plot_times)
    for s = 1:length(stripe_index)
        f_activity_mat(t,s) = mean(tr_fluo_vec(nc_stripe_id_vec==stripe_index(s)&...
            time_vec==plot_times(t)));
        v_activity_mat(t,s) =  nanmean(v_state_vec(nc_stripe_id_vec==stripe_index(s)&...
            time_vec==plot_times(t)));
    end
end
mean_fluo_fig = figure;
colormap(cm)
imagesc(f_activity_mat)
set(gca,'xtick',2:3:22,'xticklabels',1:7)
set(gca,'ytick',0:5:50)
title('average observed fluorescence over time')
saveas(mean_fluo_fig,[DataPath '\stripe_tracking_check.png'],'png')

% check to see if viterbi activity tracks stripe over time
v_state_fig = figure;
colormap(cm)
imagesc(v_activity_mat(10:end,:)-1)
set(gca,'xtick',2:3:22,'xticklabels',1:7)
set(gca,'ytick',0:5:50,'yticklabels',0:5:50)
title('fraction active viterbi states over time')
colorbar
saveas(v_state_fig,[DataPath '\viterbi_stripe_tracking_check.png'],'png')

ap_register_fig = figure;
cm = jet(128);
colormap(cm)
% colormap([cm(20,:);cm(40,:)])
stripe_bin_vec = longform_data(:,8)==4;%round(longform_data(:,8));
ap_registered = longform_data(:,5);
ap_raw = longform_data(:,4);
time_vec = longform_data(:,10);

hold on
% scatter(ap_raw,time_vec,40,stripeedit csvwrite_with_headers_bin_vec)
for i = 1:7
   s = scatter(ap_registered(longform_data(:,8)==i),time_vec(longform_data(:,8)==i),40,'MarkerfaceColor',cm(15*i,:));
   s.MarkerFaceAlpha = .05;
   s.MarkerEdgeAlpha = 0;
end
saveas(ap_register_fig,[DataPath '\registered_ap_check.png'],'png')

%%% Make Set-specific figures
set_stripe_fig = figure;
for i = 1:3
    for j = 1:4
        setID = 4*(i-1)+j;
        if setID > 11
            continue
        end
        subplot(3,4,setID);
        hold on
        for k = 1:7
            row_filter = longform_data(:,8)==k & longform_data(:,3) == setID;
            s = scatter(ap_registered(row_filter),time_vec(row_filter)/60,40,'MarkerfaceColor',cm(15*k,:));
            s.MarkerFaceAlpha = .05;
            s.MarkerEdgeAlpha = 0;
        end        
        title(['Set: ' num2str(setID)])
        xlim([.2 .9])
        ylim([0 50])
    end
end
set_stripe_fig.Position = [0 0 1024 1024];
saveas(set_stripe_fig,[DataPath '\set_specific_registered_ap.png'],'png')
% % s = scatter(ap_raw,time_vec,40,stripe_bin_vec,'filled');
% s.MarkerFaceAlpha = .1;
% s.MarkerEdgeAlpha = 0;
%% Save Data
csvwrite_with_headers([DataPath '\eve_data_longform_w_nuclei_K' num2str(K) '.csv'], ...
                       longform_data, header,9); 
                 
save([DataPath 'eve_data_longform_w_nuclei_K' num2str(K) '.mat'],'longform_data')                   