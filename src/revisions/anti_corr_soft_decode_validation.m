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
alpha = sim_struct(1).alpha;
inference_times = 20*60;
t_window = 20*60; % determines width of sliding window
project = 'anti_corr_02/';
% Set read path
out_suffix =  ['/' project '/w' num2str(w) '/states' num2str(K) '/']; 
out_prefix = '../../out/revisions/';
out_dir = [out_prefix out_suffix];

% Get list of inference files
file_list = dir([out_dir '*.mat']);
inference_results = struct;
rand_vec = [];
for f = 1:numel(file_list)
    load([out_dir file_list(f).name])
    fn = fieldnames(output);
    for i = 1:numel(fn)
        inference_results(f).(fn{i}) = output.(fn{i});
    end
    rand_vec = [rand_vec strcmp(output.sim_name,'random')];
end
rand_index = fliplr(unique(rand_vec));

%% First establish whether anticorrelation artifact is present in pure continuous traces
kon_vec_true = NaN(numel(sim_struct(1).promoter_states),2);
koff_vec_true = NaN(numel(sim_struct(1).promoter_states),2);
for rnd = 1:2
    jump_times = sim_struct(rnd).jump_times;
    promoter_states = sim_struct(rnd).promoter_states;    
    for i = 1:numel(promoter_states)
        jt = jump_times{i};
        ps = promoter_states{i};
        ps = ps(1:end-1);
        dTimes = diff(jt);
        tau_off = mean(dTimes(ps==1));
        kon_vec_true(i,rnd) = 1/tau_off/2;
        % have to be a bit more subtle for koff        
        z2 = ps==2;
        z2_diff = [diff(ps) 0];    
        out_switches = z2_diff(z2);        
        tau_2 = mean(dTimes(z2));        
        k_eff = 1/tau_2;
        gamma = sum(out_switches==-1)/sum(out_switches==1);
        koff_vec_true(i,rnd) = gamma/(gamma+1) * k_eff;
    end
%     cm = jet(128);
%     true_fig = figure;
%     scatter(kon_vec_true, koff_vec_true, 20, 'MarkerFaceColor',cm(30,:),'MarkerFaceAlpha',.2)
%     title('True Bursting Rates')
%     xlabel('k_{on}')
%     ylabel('k_{off}')
%     axis([0 .03 0 .03])
%     saveas(true_fig, [FigPath 'true_rates.png'])
end
    %% Examine Soft Decoded Array for Bias
    close all
    for i = 1:numel(rand_index)    
        inf_struct = inference_results(rand_vec==rand_index(i));
        for j = 1:numel(inf_struct)
            soft_struct = inf_struct(j).soft_struct;
            p_zz_log_soft = soft_struct.p_zz_log_soft;
            for k = 1:numel(p_zz_log_soft)
                zz = exp(p_zz_log_soft{k});
                A = nanmean(zz,3);
                A = A ./ sum(A);
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

        koff_vec_est = reshape(tr_rate_array(1,2,:),[],1);
        koff_vec_est = koff_vec_est(imag(koff_vec_est)==0); % only use non-complex rates 
        kon_vec_est = reshape(tr_rate_array(2,1,:),[],1)/2;
        kon_vec_est = kon_vec_est(imag(kon_vec_est)==0); % only use non-complex rates 
        
        [kont, si] = sort(kon_vec_true(:,i));
        kofft = koff_vec_true(si,i);
        ft = ~isnan(kofft)&~isnan(kont);
        kont = kont(ft);
        kofft = kofft(ft);
        
        p_true = polyfit(kont,kofft,2);
        p_est = polyfit(kon_vec_est,koff_vec_est,2);
        
        tr_rate_scatter = figure('Position',[0 0 1024 512]);
        subplot(1,2,1)        
        hold on
        scatter(sim_struct(i).k_on_vec,sim_struct(i).k_off_vec,10,'MarkerFaceColor',[0 0 0],...
                'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0)
        scatter(kont, kofft, 20, 'MarkerFaceColor',cm(120,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)                
        plot(kont,polyval(p_true,kont),'Color',cm(120,:),'LineWidth',2)
        
        legend('input parameters', 'induced rates')
        ylabel('k_{off}')
        xlabel('k_{on}')        
        title('Single Trace Transition Rates (Idealized Decoding)')
        axis([0 .035 0 .035])                
        
        subplot(1,2,2)
        hold on
        scatter(sim_struct(i).k_on_vec,sim_struct(i).k_off_vec,10,'MarkerFaceColor',[0 0 0],...
                'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',0)
        scatter(kon_vec_est, koff_vec_est, 20, 'MarkerFaceColor', cm(30,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
        plot(sort(kon_vec_est),polyval(p_est,sort(kon_vec_est)),'Color',cm(30,:),'LineWidth',2)        
        xlabel('k_{on}')        
        title('Single Trace Transition Rates (Actual Decoding)')
        legend('input parameters', 'induced rates')
        axis([0 .035 0 .035])        
        saveas(tr_rate_scatter,[FigPath 'tr_rate_scatter_' sim_struct(i).ID '.png'])        
    end

