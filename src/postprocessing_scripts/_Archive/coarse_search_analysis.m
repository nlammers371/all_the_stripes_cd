addpath('D:\Data\Nick\projects\all_the_stripes_cd\src\utilities\');
% Script to establish distances between sets in likelihood space

% ----------------Load Inference Results and Trace Data------------------%

% path to results
resultsPath = '../../out/eve7stripes_inf_2017_09_25/2017-09-25/mem7_states2/coarse_search/coarse_logL_search_K2G50.mat';
load(resultsPath)
% set write path
writePath = 'D:\Data\Nick\projects\all_the_stripes_cd\fig\experimental_system\eve7stripes_inf_2017_09_25\2017-09-25\mem7_states2\coarse_search\';
mkdir(writePath);

%% --------------Calculate Pairwise Distance Matix------------------------%
N = length(coarse_search_struct);
dist_mat = zeros(N,N);
for i = 1:N
    logL_i = exp(coarse_search_struct(i).logL_vec) / sum(exp(coarse_search_struct(i).logL_vec));
    coarse_search_struct(i).L_vec_norm = logL_i;
    set_i = coarse_search_struct(i).setID;
    stripe_i = coarse_search_struct(i).stripe_id;
    for j = 1:N
        logL_j = exp(coarse_search_struct(j).logL_vec) / sum(exp(coarse_search_struct(j).logL_vec));
        set_j = coarse_search_struct(j).setID;
        stripe_j = coarse_search_struct(j).stripe_id;
        c = corrcoef(logL_i,logL_j);
        dist_mat(i,j) = c(2,1);
%         dist_mat(i,j) = sum(abs(logL_i-logL_j));
    end
end
%reorder things by stripe
index_vec = 1:N;
stripe_id_vec = [coarse_search_struct.stripe_id];
[stripe_id_sorted, ranked_stripe] = sort(stripe_id_vec);
set_id_vec = [coarse_search_struct.setID];
set_id_sorted = set_id_vec(ranked_stripe);
sort_index = index_vec(ranked_stripe);
dist_mat_sorted = dist_mat(sort_index,sort_index);
% make list of tick labels
tick_list = {};
for i = 1:N
    tick_list = [tick_list{:} {[num2str(stripe_id_sorted(i)) '.' num2str(set_id_sorted(i))]}];
end
logL_space_fig = figure('Position',[0 0 1024 1024]);
colormap('winter')
imagesc(dist_mat_sorted)
set(gca,'xtick',1:30,'xticklabel', tick_list)
xtickangle(90)
set(gca,'ytick',1:30,'yticklabel', tick_list)
title('Pairwise Cross-Correlation Matrix of Parameter-Space LogL Scores')
xlabel('Stripe.Set')
saveas(logL_space_fig,[writePath 'logL_space_correlation_mat.png'],'png')
%% ------------------- Generate LongForm CSV --------------------------- %%
coarse_sorted = coarse_search_struct(ranked_stripe);
granularity = coarse_search_struct(1).granularity;
% Generate indexing array
index_array = repmat(combvec(1:granularity,1:granularity,1:granularity)',length(coarse_sorted),1);
% make longform ID vectors
stripe_vec_long = repelem([coarse_sorted.stripe_id],(granularity^3))';
set_vec_long = repelem([coarse_sorted.setID],(granularity^3))';
% longform static param arrays
mem_vec_long = repelem([coarse_sorted.w],(granularity^3))';
alpha_vec_long = repelem([coarse_sorted.alpha],(granularity^3))';
noise_vec_long = repelem([coarse_sorted.mean_noise],(granularity^3))';
% concatenate logL vectors
logL_vec_long = [coarse_sorted.logL_vec]';
L_norm_vec_long = [coarse_sorted.L_vec_norm]';
param_array = coarse_sorted(1).param_array';
param_array_rates = zeros(size(param_array));
param_array_rates(:,1) = param_array(:,1);
for i = 1:size(param_array_rates,1)
    A = zeros(2,2);
    A(2,1) = param_array(i,3);
    A(1,2) = param_array(i,2);
    A(eye(2)==1) = 1 - sum(A);
    R = logm(A)/20;
    param_array_rates(i,3) = R(2,1);
    param_array_rates(i,2) = R(1,2);
end
param_array_rates_long = repmat(param_array_rates,N,1);

%%% Combine
output = [set_vec_long stripe_vec_long param_array_rates_long index_array logL_vec_long L_norm_vec_long];
% make header
header = {'SetID', 'StripeID', 'r', 'k_on', 'k_off', 'r_Dim','k_on_Dim','k_off_Dim', 'LogLikelihood','NormLikelihood'};
%write
csvwrite_with_headers([writePath 'longform_likelihood_granularity' ...
    num2str(granularity) '.csv'], output, header);