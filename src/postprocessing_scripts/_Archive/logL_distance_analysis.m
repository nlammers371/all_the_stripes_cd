addpath('D:\Data\Nick\projects\all_the_stripes_cd\src\utilities\');
% Script to establish distances between sets in likelihood space

%-----------------Load Inference Results and Trace Data-------------------%

% set path to folder with bootstrap results
resultsFolder = 'D:\Data\Nick\projects\all_the_stripes_cd\out\eve7stripes_inf_2017_09_25\2017-09-25\inference_w7set_boots\';
% get list of files in folder
bootFiles = dir([resultsFolder 'eveSet_w7*']);
% set write path
writePath = 'D:\Data\Nick\projects\all_the_stripes_cd\fig\experimental_system\eve7stripes_inf_2017_09_25\2017-09-25\mem7_states2_set_bootstrapping\';
mkdir(writePath);
% set path to trace data
tracePath = 'D:\Data\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2017_09_25\inference_traces_t20_eve7stripes_inf_2017_09_25.mat';
% load trace data (saved as "interp_struct")
load(tracePath);

% get lists of unique set-stripe combinations
stripe_list = [];
set_list = [];
for s = 1:7
    stripe_sets = unique([interp_struct([interp_struct.stripe_id]==s).setID]);
    set_list = [set_list stripe_sets];
    stripe_list = [stripe_list linspace(s,s,length(stripe_sets))];
end
%%
setID_list = unique([interp_struct.setID]);
set_key_struct = struct;
for s = 1:length(setID_list)
    set_key_struct(s).setID = setID_list(s);
    set_struct = interp_struct([interp_struct.setID]==setID_list(s));
    set_key_struct(s).set_name = set_struct(1).set;
end
    
%%
% load inference results
inf_results = [];
for i = 1:length(bootFiles)
    load([resultsFolder bootFiles(i).name])
    inf_results = [inf_results outputs];
end
rm_sets = [inf_results.removed_set];
stripe_group = [inf_results.stripe_id];
%%
% create array to store likelihoods
K = 2;
logL_mat = zeros(length(inf_results),length(inf_results));
% outer loop determines inf set identity
for i = 1:length(inf_results)
    inf_set = inf_results(i);
    A_log = inf_set.A_log;
    v = inf_set.v;
    noise = inf_set.noise;
    pi0_log = inf_set.pi0_log;
    alpha = inf_set.alpha;
    w = inf_set.w;
    % inner iterates through all set-stripe combos
    parfor j = 1:length(inf_results)
        test_set = set_list(j);
        test_stripe = stripe_list(j);        
        %get test data
        test_data = interp_struct(([interp_struct.setID]==test_set)&...
            ([interp_struct.stripe_id]==test_stripe));
        n_dp = length([test_data.fluo]);
        fluo_values = cell(1,length(test_data));
        for k = 1:length(test_data)
            fluo_values{k} = test_data(k).fluo;
        end
        logL_mat(i,j) = likelihood_reduced_memory (fluo_values, v', noise, ...
                         pi0_log, A_log, K, w, alpha)/n_dp;
    end
end
%% Calculate Pairwise Covariance Matrix
logL_cov_mat = zeros(size(logL_mat));
logL_cells = cell(1,7);
stripe_cells = cell(1,7);
for i = 1:size(logL_mat,2)
    logL_vec_i = logL_mat(:,i);
    stripe_i = stripe_group(i);
    % get id of set removed
    rm_set_i = inf_results(i).removed_set;    
    % get list of sets used for inference (will need to be removed)
    inf_sets_i = inf_results(i).kept_sets;
    filter_vec_i = 1==(~ismember(stripe_group,stripe_i) + ...
            (stripe_group==stripe_i).*(rm_sets==rm_set_i));
    l_vec = logL_cells{stripe_i};
    logL_cells{stripe_i} = [l_vec logL_vec_i(filter_vec_i)];
    s_vec = stripe_cells{stripe_i};
    stripe_cells{stripe_i} = [s_vec stripe_group(filter_vec_i)'];
    for j = 1:size(logL_mat,2)
        logL_vec_j = logL_mat(:,j);
        stripe_j = stripe_group(j);
        % get removed set id
        rm_set_j = inf_results(j).removed_set;
        % get kept set ids
        inf_sets_j = inf_results(j);
        % make logical filter vector
        filter_vec_j = 1==(~ismember(stripe_group,stripe_j) + ...
            (stripe_group==stripe_j).*(rm_sets==rm_set_j));
        filter_vec = 1==(filter_vec_i.*filter_vec_j);
        c = corrcoef(logL_vec_i(filter_vec),logL_vec_j(filter_vec));
        logL_cov_mat(i,j) = c(2,1);
    end    
end
%% Make Relative likelihood  plots for each stripe
for i = 1:7
    stripe_L_fig = figure;
    hold on
    colormap('jet');
    cm = colormap;
    color_inc = floor(size(cm,1)/length(unique(set_list)));
    
    sets = set_list(stripe_list==i);
    stripe_id_mat = stripe_cells{i};
    stripe_logL_mat = logL_cells{i};
    stripe_mean_mat = zeros(7,length(sets));   

    legend_list = {};
    max_l = -Inf;
    min_l = Inf;
    for j = 1:size(stripe_id_mat,2)
        set_num = sets(j);
        s = stripe_id_mat(:,j);
        l = stripe_logL_mat(:,j);    
        % Normalize
        l = l-l(s==i);
        max_l = max(max(l),max_l);
        min_l = min(min(l),min_l);
        for k = 1:7
            stripe_mean_mat(k,j) = mean(l(s==k));
        end

        sctr = scatter(s,l,MarkerSize,'MarkerFaceColor',cm(color_inc*(set_num-1)+1,:),...
            'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
        set(get(get(sctr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend_list = [legend_list {['Set ' num2str(sets(j))]}];
    end
    for st = 1:length(sets)
        set_num = sets(st);
        p = plot(1:7,stripe_mean_mat(:,st),'black','LineWidth',1);
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        scatter(1:7,stripe_mean_mat(:,st),MarkerSize,...
            'MarkerFaceColor',cm(color_inc*(set_num-1)+1,:)...
            ,'MarkerEdgeColor','black')    
    end
    legend(legend_list{:},'Location','southeast')
    grid on
    title(['Normalized Log Likelihood Scores for Stripe ' num2str(i)])
    xlabel('Stripe')
    axis([0 8 min_l-.05 .05+max_l])
    
    saveas(stripe_L_fig,[writePath 'Relative_Likelihoods_Stripe' num2str(i) '.png'],'png')
    saveas(stripe_L_fig,[writePath 'Relative_Likelihoods_Stripe' num2str(i) '.eps'],'epsc')
end

%% Visualize Inference Parameters

r_2 = [inf_results.r];

rates_2 = zeros(2,length(set_list));
occupancy_2 = zeros(2,length(set_list));
for i = 1:length(set_list)
    [~, ranked_r] = sort(r_2(:,i));
    A = inf_results(i).A_mat(ranked_r,ranked_r);
    R = prob_to_rate(A,20);
    rates_2(:,i) = [R(1,2),R(2,1)];
    occupancy_2(:,i) = rates_2(:,i)./repmat(sum(rates_2(:,i)),2,1);
end

mean_rates_2 = zeros(2,7);
std_rates_2 = zeros(2,7);
mean_r_2 = zeros(2,7);
std_r_2 = zeros(2,7);
mean_occ_2 = zeros(2,7);
std_occ_2 = zeros(2,7);
for s = 1:7
   mean_rates_2(:,s) = mean(rates_2(:,stripe_group==s),2);
   std_rates_2(:,s) = std(rates_2(:,stripe_group==s),[],2);
   
   mean_r_2(:,s) = mean(r_2(:,stripe_group==s),2);
   std_r_2(:,s) = std(r_2(:,stripe_group==s),[],2);
   
   mean_occ_2(:,s) = mean(occupancy_2(:,stripe_group==s),2);
   std_occ_2(:,s) = std(occupancy_2(:,stripe_group==s),[],2);
end
MarkerSize = 75;
transition_fig = figure;
colormap('winter')
cm = colormap;
hold on
for k = 1:K
    p = plot(1:7,mean_rates_2(k,:),'black','LineWidth',1.5);
end
sctr = scatter(repelem(stripe_group,2), reshape(rates_2,1,[]),MarkerSize, ...
    cm(repmat([15,45],1,length(set_list)),:),'filled');
set(sctr,'MarkerFaceAlpha',.3)
sctr = scatter(repelem(1:7,2), reshape(mean_rates_2,1,[]),MarkerSize, ...
    cm(repmat([15,45],1,7),:),'filled');
set(sctr,'MarkerFaceAlpha',.5,'MarkerEdgeColor','black')
grid on
ylabel('s^{-1}')
xlabel('Stripe')
title('Transition Rates by Stripe')
saveas(transition_fig,[writePath 'Transition_Rates.png'],'png')
saveas(transition_fig,[writePath 'Transition_Rates.eps'],'epsc')
%%%Loading Rates
r_fig = figure;
hold on
for k = 1:K
    p = plot(1:7,mean_r_2(k,:),'black','LineWidth',1.5);
end
sctr = scatter(repelem(stripe_group,2), reshape([inf_results.r],1,[]),MarkerSize, ...
    cm(repmat([15,45],1,length(set_list)),:),'filled');
set(sctr,'MarkerFaceAlpha',.3)
sctr = scatter(repelem(1:7,2), reshape(mean_r_2,1,[]),MarkerSize, ...
    cm(repmat([15,45],1,7),:),'filled');
set(sctr,'MarkerFaceAlpha',.5,'MarkerEdgeColor','black')
grid on
ylabel('AU s^{-1}')
xlabel('Stripe')
title('Loading Rates by Stripe')
saveas(r_fig,[writePath 'Loading_Rates.png'],'png')
saveas(r_fig,[writePath 'Loading_Rates.eps'],'epsc')
%%%occupancy
occ_fig = figure;
hold on
for k = 1:K
    p = plot(1:7,mean_occ_2(k,:),'black','LineWidth',1.5);
end
sctr = scatter(repelem(stripe_group,2), reshape(occupancy_2,1,[]),MarkerSize, ...
    cm(repmat([15,45],1,length(set_list)),:),'filled');
set(sctr,'MarkerFaceAlpha',.3)
sctr = scatter(repelem(1:7,2), reshape(mean_occ_2,1,[]),MarkerSize, ...
    cm(repmat([15,45],1,7),:),'filled');
set(sctr,'MarkerFaceAlpha',.5,'MarkerEdgeColor','black')
grid on
ylabel('Share')
xlabel('Stripe')
title('Occupancy Shares by Stripe')
saveas(occ_fig,[writePath 'Occupancy.png'],'png')
saveas(occ_fig,[writePath 'Occupancy.eps'],'epsc')

%% Visualize Covariance Matrix
cov_fig = figure('Position',[0 0 1024 1024]);
colormap('winter')
imagesc(logL_cov_mat)
set(gca,'xtick',1:39,'xticklabel', stripe_group)
set(gca,'ytick',1:39,'yticklabel', stripe_group)
h = colorbar;
title('Pairwise Correlation in Log Likelihood Scores')
saveas(cov_fig,[writePath 'LogL_Correlation_Matrix.png'],'png')
saveas(cov_fig,[writePath 'LogL_Correlation_Matrix.eps'],'epsc')
