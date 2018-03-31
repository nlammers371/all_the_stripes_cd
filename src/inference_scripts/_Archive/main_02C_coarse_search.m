addpath('../utilities')
%Script to perform coarse-grained exploration of Likelihood Space for a
%specified architecture
% -------------------- Read in Inference Files -------------------------%%
%memory assumed for inference
w = 7;
%state count 
K = 2;
%Tres
dT = 20;
%set bin range
region_range = [1:7];
% ID variables
date_str = '2017-09-25';
project = 'eve7stripes_inf_2017_09_25';

folder_path = ['../../out/' project '/' date_str '/inference_w' num2str(w) '/individual_results/'];
files = dir(folder_path);
% collect filenames fitting specified inf type params
filenames = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
            ~isempty(strfind(files(i).name,['K' num2str(K)])) 
        filenames = [filenames {files(i).name}];
    end    
end

%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 1;
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'    
    try
        load([folder_path filenames{f}], 'output');
        for fn = fieldnames(output)'
            glb_all(f_pass).(fn{1}) = output.(fn{1});
        end
        f_pass = f_pass + 1;
    catch
        warning(['File ' num2str(f) ' failed to load'])
    end
end
% filter results
glb_all = glb_all(ismember([glb_all.bin],region_range));
% set writepath
OutPath = ['../../out/' project '/' date_str '/mem' ...
    num2str(w) '_states' num2str(K) '/coarse_search/' ];
mkdir(OutPath);

%load traces (saved as "interp_struct")
load(['../../dat/' project '/inference_traces_t' num2str(dT) '_' project '.mat']);
interp_struct = interp_struct(ismember(floor([interp_struct.stripe_id]),region_range));

%% -------------------Find Bounds of Parameter Space ------------------- %%
% initiation rates
r_mat = zeros(K,2);
r_mat(:,1) = Inf;
r_mat(:,2) = -Inf;
% transition probabilities
A_mat = zeros(K,K,2);
A_mat(:,:,1) = Inf;
A_mat(:,:,2) = -Inf;
% noise
noise_vec = [Inf,-Inf];
% memory
mem_vec = w;
% iterate through sets
for i = 1:length(glb_all)
    % loading rates
    [r, ranked_r] = sort(glb_all(i).r);
    r_mat(:,1) = min([r_mat(:,1),r],[],2);
    r_mat(:,2) = max([r_mat(:,2),r],[],2);
    % transition probs
    A = glb_all(i).A_mat;
    A = A(ranked_r,ranked_r);
    A_mat(:,:,1) = min(cat(3, A, A_mat(:,:,1)),[],3);
    A_mat(:,:,2) = max(cat(3, A, A_mat(:,:,2)),[],3);
    % noise
    noise = glb_all(i).noise;
    noise_vec(1) = min([noise_vec(1),noise]);
    noise_vec(2) = min([noise_vec(2),noise]);
end
%% --------------------Divide Traces by Set and Stripe------------------ %%
unique_sets_vec = unique([interp_struct.setID]);
all_sets_vec = [];
all_stripes_vec = [];
min_dp = 500;
for s = 1:length(unique_sets_vec)
    set_stripes = unique([interp_struct([interp_struct.setID]==unique_sets_vec(s)).stripe_id]);
    %filter out sets with insufficient data
    inf_stripes = [];
    inf_sets = [];
    for st = set_stripes
        N = length([interp_struct(([interp_struct.setID]==unique_sets_vec(s))&...
                    ([interp_struct.stripe_id]==st)).time]);
        if N >= min_dp
            inf_stripes = [inf_stripes  st];
            inf_sets = [inf_sets unique_sets_vec(s)];
        end
    end
    all_stripes_vec = [all_stripes_vec inf_stripes];
    all_sets_vec = [all_sets_vec inf_sets];
end

%% --------------------- Generate Search Space --------------------------%%
% set number of marks along each axis
granularity = 10;
% fix noise for now
mean_noise = mean([glb_all.noise]);
% fix pi0 PDF
mean_pi0 = mean(reshape([glb_all.pi0],K,[]),2);
%set alpha
alpha = glb_all(1).alpha;

% define parameter vectors
param_cell = cell(1,K^2-1);
for r_ind = 1:K-1
    param_cell{r_ind} = linspace(r_mat(r_ind+1,1),r_mat(r_ind+1,2),granularity);
end
%track indices
rc_mat = zeros(K^2-K,2);
s_vec = 1:K;
for A_ind = K:K^2-1
    i = A_ind - K + 1;    
    row = ceil(i/(K-1));
    ind_vec = s_vec(s_vec~=row);
    col = ind_vec(i-(K-1)*(row-1));
    rc_mat(i,:) = [row, col];    
    param_cell{A_ind} = linspace(A_mat(row,col,1),A_mat(row,col,2),granularity);
end
%% ------------------Conduct Search of Parameter Space------------------ %%
coarse_search_struct = struct;
%generate master param array
master_p_array = combvec(param_cell{:});
for i = 1:length(all_sets_vec)
    inf_struct = interp_struct(([interp_struct.setID]==all_sets_vec(i))&...
                ([interp_struct.stripe_id]==all_stripes_vec(i))); 
    ndp = length([inf_struct.fluo]);
    fluo_values = cell(1,length(inf_struct));
    for tr = 1:length(inf_struct)
        fluo_values{tr} = inf_struct(tr).fluo;
    end
    coarse_search_struct(i).setID = all_sets_vec(i);
    coarse_search_struct(i).stripe_id = all_stripes_vec(i);
    coarse_search_struct(i).ndp = ndp;
    coarse_search_struct(i).granularity = granularity;
    coarse_search_struct(i).K = K;
    coarse_search_struct(i).w = w;
    coarse_search_struct(i).alpha = alpha;
    coarse_search_struct(i).mean_pi0 = mean_pi0;
    coarse_search_struct(i).mean_noise = mean_noise;
    % This is really inefficient
    coarse_search_struct(i).param_array = master_p_array;
    % Make array to store results
    logL_vec = zeros(1,size(master_p_array,2));
    % Calculate likelihoods    
    for j = 1:size(master_p_array,2)
        param_vec = master_p_array(:,j);
        r_inf = zeros(K,1);
        for r_ind = 1:K-1
            r_inf(r_ind+1) = param_vec(r_ind);
        end
        v_inf = r_inf*dT;
        A_inf = zeros(K,K);
        for A_ind = K:K^2-1
            ind = A_ind - K + 1;
            A_inf(rc_mat(ind,1),rc_mat(ind,2)) = param_vec(A_ind);
        end
        A_inf(eye(K)==1) = 1 - sum(A_inf);        
        A_log = log(A_inf);
        logL_normed = likelihood_reduced_memory (fluo_values, v_inf', mean_noise, ...
                         log(mean_pi0), A_log, K, w, alpha)/ndp;
        logL_vec(j) = logL_normed;
    end    
    coarse_search_struct(i).logL_vec = logL_vec;
end
save([OutPath 'coarse_logL_search_K' num2str(K) 'G' num2str(granularity) '.mat'],'coarse_search_struct')


