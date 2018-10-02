% script to generate data sets for testing possible link between detection
% threshold and apparent transition rate correlation
clear 
close all
% set paths
addpath('../../../hmmm/src/utilities');
ResultsPath = '..\..\dat\_initial_eLife_submission\eve7stripes_inf_2018_04_28\w7_t20_alpha14_f1_cl1_no_ends1_tbins1\K3_summary_stats\';
DataPath = '../../dat/revisions/anti_corr/';
mkdir(DataPath);
load([ResultsPath 'hmm_results_t_window30_t_inf40.mat'])
FigPath = '../../fig/revisions/anti_corr/validation_03/';
mkdir(FigPath);
% extract average bursting parameters
bin_index = [hmm_results.binID];
% set simulation parameteres
K = 3; % N states
w = 5; % memory
r_mean = hmm_results(bin_index==0).initiation_mean/60; % initiation rates
r_mean = [0 r_mean(2) 2*r_mean(2)];
pi0 = [1 1 1]/3; % initial state PDF
alpha = w/7*hmm_results(bin_index==0).alpha; %MS2 rise time
dT = hmm_results(bin_index==0).dT; % time resolution
noise_mean = w/7*hmm_results(1).noise_mean; % noise level (AU)
seq_length = 60/dT*40;
trace_time = (1:seq_length)*20;

% set simulation parameters
n_traces = 1000; % number of traces per simulation

constant_rates = [.01 .02 .03]; % constant rate magnitudes to test
sim_rate_bounds = [.005 .05]; % range of rate sto simulate

sim_ids = {'koff_constant','kon_constant'};
sim_struct = struct;
% generate sample traces
for i = 1:numel(sim_ids)
    for j = 1:numel(constant_rates)
        ID = sim_ids{i};
        ID_iter = ([ID '_' num2str(constant_rates(j))]);
        c_rate = constant_rates(j);
        
        ind = (i-1)*numel(constant_rates) + j;
        
        traces = cell(n_traces,1);
        trajectories = cell(n_traces,1);
        jumps = cell(n_traces,1);
        f_vec = NaN(1,n_traces);        
        if strcmp(ID,'koff_constant')
            koff_vec = repelem(c_rate,n_traces);
            kon_vec = rand(size(koff_vec))*(sim_rate_bounds(2)-sim_rate_bounds(1))+sim_rate_bounds(1);
        elseif strcmp(ID,'kon_constant')
            kon_vec = repelem(c_rate,n_traces);
            koff_vec = rand(size(koff_vec))*(sim_rate_bounds(2)-sim_rate_bounds(1))+sim_rate_bounds(1);
        end
        for n = 1:n_traces                 
            kon = kon_vec(n);
            koff = koff_vec(n);
            % make rate matrix
            R = [-2*kon koff 0 ; 2*kon -koff-kon 2*koff; 0 kon -2*koff];            
            gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                            K, w, R, dT, r_mean, noise_mean, pi0);
            traces{n} = gillespie.fluo_MS2;
            if sum(isnan(gillespie.fluo_MS2))>0
                error('afa')
            end
            trajectories{n} = gillespie.naive_states;
            jumps{n} = gillespie.transition_times;
            kon_vec(n) = kon;
            koff_vec(n) = koff;
            f_vec(n) = nanmean(traces{n});        
        end
        sim_struct(ind).ID = sim_ids{i};    
        sim_struct(ind).dT = dT;
        sim_struct(ind).fluo_data  = traces;
        sim_struct(ind).trace_time = trace_time;
        sim_struct(ind).promoter_states  = trajectories;
        sim_struct(ind).jump_times = jumps;    
        sim_struct(ind).r_emission = r_mean;
        sim_struct(ind).noise = noise_mean;
        sim_struct(ind).w = w;
        sim_struct(ind).alpha = alpha;
        sim_struct(ind).K = K;
        sim_struct(ind).rate_bounds = sim_rate_bounds;               
        sim_struct(ind).k_off_vec = koff_vec;
        sim_struct(ind).k_on_vec = kon_vec;

        rate_scatter = figure; 
        if strcmp(ID,'koff_constant')
            [~, si ] = sort(kon_vec);
        else
            [~, si ] = sort(koff_vec);
        end
        scatter(kon_vec(si),koff_vec(si),40,f_vec,'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
        xlim(sim_rate_bounds)
        ylim(sim_rate_bounds)
        xlabel('k_{on} (s^{-1})')
        ylabel('k_{off} (s^{-1})')
        title(['k_{on} vs. k_{off} (' sim_ids{i} ')'])
        colorbar
        saveas(rate_scatter,[FigPath 'rate_scatter_' ID_iter '.png'])
    end
end
save([DataPath 'sim_data_03.mat'],'sim_struct')