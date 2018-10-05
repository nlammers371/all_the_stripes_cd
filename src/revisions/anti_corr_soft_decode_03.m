% Conduct analyses to assess impact of detection limits on apparent
% anticorrelation of kon and koff
close all
clear 
addpath('../../../hmmm/src/utilities');

% Load simulation data
DataPath = '../../dat/revisions/anti_corr/';
load([DataPath 'sim_data_03.mat'])
% Figure write path
FigPath = '../../fig/revisions/anti_corr/cycle_time/';
mkdir(FigPath);
% Load inference results
K = sim_struct(1).K; % N States
w = sim_struct(1).w; % Memory
dT = sim_struct(1).dT;
alpha = sim_struct(1).alpha;
inference_times = 20*60;
t_window = 20*60; % determines width of sliding window
project = 'anti_corr_03/';
% Set read path
out_suffix =  ['/' project '/w' num2str(w) '/states' num2str(K) '/']; 
out_prefix = '../../out/revisions/';
out_dir = [out_prefix out_suffix];

% Get list of inference files
file_list = dir([out_dir '*.mat']);
inference_results = struct;
ID_cell = {};
inf_id_vec = [];
for f = 1:numel(file_list)
    load([out_dir file_list(f).name])
    fn = fieldnames(output);
    for i = 1:numel(fn)
        inference_results(f).(fn{i}) = output.(fn{i});
    end    
    ID_cell = [ID_cell{:} {output.sim_name}];
    inf_id_ind = strfind(file_list(f).name,'_id');
    inf_id = str2double(file_list(f).name(inf_id_ind+3:end-4));
    inf_id_vec = [inf_id_vec inf_id];
end
[inf_id_vec, si] = sort(inf_id_vec);
inference_results = inference_results(si);
sim_index = 1:numel(sim_struct);
inf_sim_id = repelem(sim_index,5);

%% Examine Soft Decoded Array for Bias
% try average two different ways: 1) decode traces individually and average
% 2) average trace tranisition probs by "region" then decode
smothing_window = .00025; % assuming we have proxy variable for ref var
close all
boot_size = numel(inference_results(1).boot_indices);
results_struct = struct;
off_table = [];
on_table = [];

sim_traces = [];
sim_time = [];
sim_kon = [];
sim_koff = [];
inferred_kon = [];
inferred_koff = [];
for i = 1:numel(sim_index) 
    ID = sim_struct(i).ID;
    inf_struct = inference_results(inf_sim_id==sim_index(i));    
    s_struct = sim_struct(i);
    if strfind(ID,'koff')
        ref_var = 'k_on_vec';                
        response_var = 'k_off_vec';                
        indices = [1,2];
        factor = 1;
    else
        ref_var = 'k_off_vec';
        response_var = 'k_on_vec';                
        indices = [2,1];
        factor = 2;
    end    
    ref_vec = s_struct.(ref_var);    
    A_array = NaN(K,K,numel(inf_struct)*boot_size);
    ref_var_vec = NaN(1,numel(inf_struct)*boot_size);
    boot_indices_full = [];    
    for j = 1:numel(inf_struct)
        % record true ref rates
        bi = inf_struct(j).boot_indices;
        ref_var_vec((j-1)*boot_size+1:j*boot_size) = ref_vec(bi);
        boot_indices_full = [boot_indices_full bi];        
        sim_traces = [sim_traces [sim_struct(i).fluo_data{bi}]];
        sim_time = [sim_time [sim_struct(i).fluo_data{bi}]];
        sim_kon = [sim_kon [sim_struct(i).k_on_vec(bi)]];
        sim_koff = [sim_koff [sim_struct(i).k_off_vec(bi)]];
        
        % iterate through decoded traces
        soft_struct = inf_struct(j).soft_struct;
        p_zz_log_soft = soft_struct.p_zz_log_soft;        
        for k = 1:numel(p_zz_log_soft)
            zz = exp(p_zz_log_soft{k});
            A = nanmean(zz,3);
            A = A ./ sum(A);
            A_array(:,:,(j-1)*boot_size+k) = A;            
            if k == 1 && j == 1
                tr_prob_array = A;
                tr_rate_array = logm(A)/dT;
            else
                tr_prob_array = cat(3,tr_prob_array, A);
                tr_rate_array = cat(3,tr_rate_array,logm(A)/dT);
            end
        end
    end    
    cm = jet(128);
    [ref_var_sorted, si] = sort(ref_var_vec+rand(size(ref_var_vec))/1e6);
    kon_vec_inf = reshape(tr_rate_array(2,1,:),1,[])/2;
    koff_vec_inf = reshape(tr_rate_array(1,2,:),1,[]);
    
    im_ft = imag(kon_vec_inf)==0&imag(koff_vec_inf)==0;
    kon_vec_inf(~im_ft) = NaN;
    koff_vec_inf(~im_ft) = NaN;
    kon_vec_inf = real(kon_vec_inf);
    koff_vec_inf = real(koff_vec_inf);
    % record
    inferred_kon = [inferred_kon kon_vec_inf];
    inferred_koff = [inferred_koff koff_vec_inf];
    
    tr_rate_vec = reshape(tr_rate_array(indices(1),indices(2),si),1,[])/factor; % sort rate array
    im_vec = imag(tr_rate_vec)>0;
    tr_rate_vec = real(tr_rate_vec);
    tr_rate_vec(im_vec) = NaN;
    A_array = A_array(:,:,si); % sort A array
    freq = round(1/(diff(sim_struct(i).rate_bounds)/size(A_array,3)));
    % take average response as a function of true underlying ref var     
    [tr_rates_resamp, r_ref_resamp] = resample(tr_rate_vec,ref_var_sorted,freq);        
    tr_rates_smooth = imgaussfilt(tr_rates_resamp,25);    
    for m = 1:K
        for n = 1:K
            [a_vec, a_ref_resamp] = resample(reshape(A_array(n,m,:),1,[]),ref_var_sorted,freq);
            if n==1 && m==1
                A_resamp = NaN(K,K,numel(a_ref_resamp));
                A_smooth = NaN(size(A_resamp));                
            end
            A_resamp(n,m,:) = a_vec;            
            A_smooth(n,m,:) = imgaussfilt(a_vec,25);
        end
    end
    % calculate rates from A arrays
    tr_rates_agg_raw = NaN(1,numel(a_ref_resamp));
    tr_rates_agg_smooth = NaN(1,numel(a_ref_resamp));
    for m = 1:numel(tr_rates_agg_raw)
        R_raw = logm(A_resamp(:,:,m))/dT;
        tr_rates_agg_raw(m) = R_raw(indices(1),indices(2))/factor;
        R_smooth = logm(A_smooth(:,:,m))/dT;
        tr_rates_agg_smooth(m) = R_smooth(indices(1),indices(2))/factor;
    end
    true_response = unique(sim_struct(i).(response_var));
    cycle_rate = (r_ref_resamp.*true_response) ./ (r_ref_resamp+true_response);
    cycle_rate_raw = (ref_var_sorted*true_response) ./ (ref_var_sorted+true_response);
    tr_rate_scatter = figure('Position',[0 0 1024 512]);             
    hold on
    scatter(cycle_rate_raw,tr_rate_vec,10,'MarkerFaceColor',[0 0 0],...
            'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',0)
    plot(cycle_rate,tr_rates_smooth,'Color',cm(25,:),'LineWidth',1.5);
    plot(cycle_rate,tr_rates_agg_smooth,'Color',cm(110,:),'LineWidth',1.5);    
    plot(cycle_rate,repelem(true_response,numel(a_ref_resamp)),'--','Color','black')
    legend('single trace estimates', 'rate averaging', 'probability averaging',...
        'ground truth','Location','southwest')
    if strfind(ID,'koff')
        ylabel('k_{off}')           
    else
        ylabel('k_{on}')                
    end    
    xlabel('cycle freq (s^{-1})')
    title('Inference Accuracy as a Function of Cycle Frequency')
    ylim([0 1.5*max(tr_rates_agg_smooth)])                
    grid on
    saveas(tr_rate_scatter,[FigPath 'tr_rate_scatter_' response_var '_' ...
        num2str(round(true_response*1000)) '.png'])        
    close all
    if strfind(ID,'koff')
        ref_var = 'k_on_true';
        resp_base = 'k_off';
    else
        ref_var = 'k_off_true';
        resp_base = 'k_on';
    end
    table = array2table([double(repelem(i,numel(a_ref_resamp))')...
        a_ref_resamp' repelem(true_response,numel(a_ref_resamp))' cycle_rate'...
        tr_rates_resamp' tr_rates_smooth' tr_rates_agg_smooth'],'VariableNames',...
        {'sim_id' ref_var [resp_base '_true'] 'cycle_freq'...
        [resp_base '_raw'] [resp_base '_avg_rates'] [resp_base '_avg_probs']});        
        
    if strfind(ID,'koff')
        off_table = [off_table ; table];        
    else
        on_table = [on_table ; table];        
    end
end

traces = array2table([repelem(1:numel(sim_struct),numel(sim_time)/numel(sim_struct))'...
    sim_time' sim_traces' repelem(sim_kon,120)', repelem(sim_koff,120)' ...
    repelem(inferred_kon,120)' repelem(inferred_koff,120)'],...
    'VariableNames', {'sim_id' 'time' 'fluo' 'k_on_true' 'k_off_true' 'k_on_inf' 'k_off_inf'});

writetable(on_table,[DataPath 'k_on_simulations.csv'])
writetable(off_table,[DataPath 'k_off_simulations.csv'])
writetable(traces,[DataPath 'sim_traces.csv'])