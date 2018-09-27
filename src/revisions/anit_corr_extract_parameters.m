% Conduct analyses to assess impact of detection limits on apparent
% anticorrelation of kon and koff
close all
clear 
addpath('../../../hmmm/src/utilities');

% Load simulation data
DataPath = '../../dat/revisions/anti_corr/';
load([DataPath 'sim_data_01.mat'])
% Figure write path
FigPath = '../../fig/revisions/anti_corr/';
mkdir(FigPath);
% Load inference results
K = sim_struct(1).K; % N States
w = sim_struct(1).w; % Memory
dT = sim_struct(1).dT;
load([DataPath '/viterbi_fits_w' num2str(w) '_K' num2str(K) '.mat']);
% Extract indexing vectors
random_vec = [viterbi_fit_struct.random];
random_index = unique(random_vec);
sigma_vec = [viterbi_fit_struct.sigma_thresh];
sigma_index = unique(sigma_vec);
% Generate scatter plots to compare inferred relations with actual
results_struct = struct;
mkdir([FigPath '/results_scatters_rnd/'])
for i = 1:numel(random_index)    
    rnd = random_index(i);
    % extract "true" bursting parameters from promoter trajectories
    promoter_states = sim_struct(i).promoter_states;
    jump_times = sim_struct(i).jump_times;
    kon_true = NaN(1,numel(jump_times));
    koff_agg_true = NaN(1,numel(jump_times));
    koff_alt_true = NaN(1,numel(jump_times));
%     trace_time = sim_struct(i).trace_time;
    trace_time = 1:1:max(sim_struct(i).trace_time);
    for n = 1:numel(jump_times)
        p_vec = promoter_states{n};
        j_vec = jump_times{n};
        p_vec_discrete = NaN(size(trace_time)); 
        for t = 1:numel(trace_time)
            t_index = find(trace_time(t)>=j_vec,1,'last');
            if isempty(t_index)
                error('sigh')
            end
            p_vec_discrete(t) = p_vec(t_index);
        end
        zb = p_vec_discrete>1; % convert to binary vec
        zb_diff = [0 diff(zb)]; % find switches
        switch_indices = find(zb_diff); 
        state_ids = zb(switch_indices);
        state_ids = state_ids(1:end-1);
        sw_dT = diff(switch_indices)*1;
        off_durations = sw_dT(state_ids==0);
        on_durations = sw_dT(state_ids==1);            
        % save
        kon_true(n) = 1/mean(off_durations)/2;
        koff_agg_true(n) = 1/mean(on_durations);            
        % now perform alternative koff calc for comparison
        z2 = p_vec_discrete==2;
        z2_diff = [diff(p_vec_discrete) 0];
        out_switches = z2_diff(z2);
        switch_indices = find([0 z2_diff(1:end-1)]);
        state_ids = z2(switch_indices);
        state_ids = state_ids(1:end-1);
        sw_dT = diff(switch_indices)*1;
        two_durations = sw_dT(state_ids==1);            
        k_eff = 1/mean(two_durations);
        gamma = sum(out_switches==-1)/sum(out_switches==1);
        koff_alt_true(n) = gamma/(gamma+1) * k_eff;
    end
    for j = 1:numel(sigma_index)
        sgm = sigma_index(j);
        ind = (i-1)*numel(sigma_index) + j;
        v_fits = viterbi_fit_struct(sigma_vec==sgm&random_vec==rnd).v_fits;
        kon_inf = NaN(1,numel(v_fits));
        koff_agg_inf = NaN(1,numel(v_fits));
        koff_alt_inf = NaN(1,numel(v_fits));
        f_vec = NaN(1,numel(v_fits));
        for k = 1:numel(v_fits)
            % estimate trace-specific bursting paramters from viterbi fits
            z_viterbi = v_fits(k).z_viterbi;
            zb = z_viterbi>1; % convert to binary vec
            zb_diff = [0 diff(zb)]; % find switches
            switch_indices = find(zb_diff); 
            state_ids = zb(switch_indices);
            state_ids = state_ids(1:end-1);
            sw_dT = diff(switch_indices)*dT;
            off_durations = sw_dT(state_ids==0);
            on_durations = sw_dT(state_ids==1);            
            % save
            kon_inf(k) = 1/mean(off_durations)/2;
            koff_agg_inf(k) = 1/mean(on_durations);
            f_vec(k) = nanmean(v_fits(k).fluo_thresh);            
            % now perform alternative koff calc for comparison
            z2 = z_viterbi==2;
            z2_diff = [diff(z_viterbi) 0];
            out_switches = z2_diff(z2);
            switch_indices = find([0 z2_diff(1:end-1)]);
            state_ids = z2(switch_indices);
            state_ids = state_ids(1:end-1);
            sw_dT = diff(switch_indices)*dT;
            two_durations = sw_dT(state_ids==1);            
            k_eff = 1/mean(two_durations);
            gamma = sum(out_switches==-1)/sum(out_switches==1);
            koff_alt_inf(k) = gamma/(gamma+1) * k_eff;
        end        
        % get actual param values
        kon_base = sim_struct(i).k_on_vec;
        koff_base = sim_struct(i).k_off_vec;
        % compare
        cm = jet(128);
        comp_fig = figure('Position',[0 0 1024 1024]);
        subplot(2,2,1)
        hold on
        scatter(kon_inf,koff_agg_inf,40,'MarkerFaceColor',cm(30,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
        scatter(kon_inf,koff_alt_inf,40,'MarkerFaceColor',cm(120,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
        legend('method 1', 'method 2')
        title(['Inferred Relationship (Threshold: ' num2str(round(sgm)) ')'])
        ylabel('k_{off}')
        xlabel('k_{on}')
        axis([0 .02 0 .02])
        
        subplot(2,2,2)
        hold on        
%         scatter(kon_true,koff_agg_true,40,'MarkerFaceColor',cm(30,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);        
        scatter(kon_true,koff_alt_true,40,'MarkerFaceColor',cm(120,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);        
%         scatter(kon_base,koff_base,40,'MarkerFaceColor',cm(60,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0);
        legend('method 1', 'method 2', 'ground truth')
        title('True Relationship')
        ylabel('k_{off}')
        xlabel('k_{on}')   
        axis([0 .04 0 .04])
        
        subplot(2,2,3)
        hold on
        scatter(kon_base,kon_true,40,'MarkerFaceColor',cm(30,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
        title(['Inferred vs. Actual k_{on} (Threshold: ' num2str(round(sgm)) ')'])
        ylabel('k_{on} (inferred)')
        xlabel('k_{on} (actual)')
        axis([0 .02 0 .02])
        
        subplot(2,2,4)
        hold on
        scatter(koff_base,koff_alt_true,40,'MarkerFaceColor',cm(30,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
%         scatter(koff_agg_true,koff_alt_inf,40,'MarkerFaceColor',cm(120,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
%         legend('method 1', 'method 2')
        title(['Inferred vs. Actual k_{off} (Threshold: ' num2str(round(sgm)) ')'])
        ylabel('k_{off} (inferred)')
        xlabel('k_{off} (actual)')
        axis([0 .02 0 .02])                
        error('afa')
        saveas(comp_fig,[FigPath '/results_scatters_rnd/' num2str(rnd) '_sgm' num2str(round(sgm)) '.png'])        
        results_struct(ind).random = rnd;
        results_struct(ind).threshold = sgm;
        results_struct(ind).kon_true = kon_true;
        results_struct(ind).koff_true = koff_agg_true;
        results_struct(ind).kon_inf = kon_inf;
        results_struct(ind).koff_inf = koff_agg_inf;
        
        close all
    end
end