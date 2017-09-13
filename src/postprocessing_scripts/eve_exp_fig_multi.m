addpath('../utilities/');

%------------------------------Set System Params--------------------------%
%memory assumed for inference
w = 8;
%states used for final inference
K = 2;
%set bin range
bin_range = [1:2];
datatype = 'weka';
date_str = '2017-09-08_test';
project = ['eve7stripes_inf'];
interp_data_name = ['_w' num2str(w)];
raw_data_name = [''];
folder_path = ['../../out/' project '/' date_str '/inference' interp_data_name '/individual_results/'];
files = dir(folder_path);
%Include original (un-fit) transition rates in plots?
plot_orig_rates = 0;
filenames = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ~isempty(strfind(files(i).name,['K' num2str(K)]))
        filenames = [filenames {files(i).name}];
    end
end

if isempty(filenames)
    error('No file with specified date string, state count, and memory found')
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

OutPath = ['../../fig/experimental_system' '/' project '/' date_str '/mem' ...
    num2str(w) '_states' num2str(K) '/' ];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

%load traces (saved as "interp_struct")
traces_all = load(['../../dat/' project '/inference_traces' interp_data_name '_' project '.mat']);
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.stripe_id]),bin_range));

dT = traces_all(1).dT;
%% Perform Rate Fitting As Necessary
%Define convenience Arrays and vectors
glb_bin_index = [glb_all.bin];
alpha = glb_all(1).alpha;
bin_index = [glb_all.bin];
bin_vec = repelem(bin_index,1,K);
bin_vec_u = unique(bin_index);
colormap('jet');
cmap = colormap;
if K == 3
    cm = flipud([cmap(50,:) ; cmap(30,:); cmap(20,:)]);
    legend_string = [{'Off'}, {'On1'},{'On2'}];
elseif K == 2
    cm = flipud([cmap(30,:); cmap(20,:)]);
    legend_string = [{'Off'}, {'On1'}];
end
c_mat = repmat(cm,length([glb_all.bin]),1);
%Perform Rate Fitting
dwell_times = zeros(K,length(glb_all));
for i = 1:length(glb_all)    
    [~, ranked_v] = sort(glb_all(i).v);
    A = reshape(glb_all(i).A,K,K);
    A = A(ranked_v, ranked_v);
    %Obtain raw R matrix
    R = prob_to_rate(A,dT);
    glb_all(i).R_mat = R(ranked_v, ranked_v);
    Rcol = reshape(R(ranked_v, ranked_v),1,[]);
    glb_all(i).R = Rcol;
    %Check for imaginary and negative elements. If present, perform rate
    %fit
    if ~isreal(Rcol)||(sum(Rcol<0)>K)
        out = prob_to_rate_fit_sym(A, dT, 'gen', .005, 1);            
        glb_all(i).R_fit = out.R_out;
        r_diag = diag(out.R_out);
    else
        glb_all(i).R_fit = glb_all(i).R_mat;
        r_diag = diag(glb_all(i).R_mat); 
    end
    dwell_times(:,i) = -r_diag.^-1;
end
%% Plot Dwell Times
avg_dwell = zeros(K, length(bin_vec_u));
std_dwell = zeros(K, length(bin_vec_u));
for i = 1:length(bin_vec_u)
    bin = bin_vec_u(i);
    for k = 1:K
        avg_dwell(k,i) = mean(dwell_times(k,bin_index==bin));        
        std_dwell(k,i) = std(dwell_times(k,bin_index==bin));        
    end
end
%Make Dwell time figure
dwell_fig = figure;
hold on
for k = 1:K
    errorbar(bin_vec_u, avg_dwell(k,:),std_dwell(k,:),'LineWidth',1,'Color','black');
end
sctr_dwell = scatter(bin_vec', reshape(dwell_times,1,[])', 75, c_mat,...
    'o', 'filled', 'MarkerFaceAlpha', .2);
sctr_emission_avg = scatter(repelem(bin_vec_u,1,K), reshape(avg_dwell,...
    1,[]) , 75, repmat(cm,length(bin_vec_u),1), 's','filled',...
    'MarkerEdgeColor', 'black','MarkerFaceAlpha',1);
% legend(sctr_emission_avg(:) ,legend_string{:})
set(gca,'FontName','Lucida Sans Regular')
axis([(min([glb_all.bin])-1) (max([glb_all.bin])+1) 0 1.1*max(max(dwell_times))])
title('Dwell Times by AP Position');
xlabel('Relative Position');
ylabel('Dwell Time (s)');
grid on 
saveas(dwell_fig, [OutPath '/dwell_times.eps'], 'epsc');
saveas(dwell_fig, [OutPath '/dwell_times.png'], 'png');

%% Generate Occupancy Plot
occupancy_ap = zeros(K,length(glb_all));
occupancy_ap_v = zeros(K,length(glb_all));
emission_rates = zeros(K,length(glb_all));
for i = 1:length(glb_all)
    [emission_rates(:,i), ranked_r] = sort([glb_all(i).r]);
    R = glb_all(i).R_fit;                                       
    [V,D] = eig(R);
    steady = V(:,1)./sum(V(:,1));
    occupancy_ap(:,i) = steady;
    glb_all(i).occupancy = steady;
end
init_vec = reshape(emission_rates,1,[]);

%Calculate averages
pi0_mat = vertcat(glb_all.pi0)';
avg_initiation = zeros(K,length(bin_vec_u));
std_initiation = zeros(K,length(bin_vec_u));
avg_pi0 = zeros(K,length(bin_vec_u));
std_pi0 = zeros(K,length(bin_vec_u));
avg_occupancy = zeros(K,length(bin_vec_u));
std_occupancy = zeros(K,length(bin_vec_u));
avg_noise = zeros(1,length(bin_vec_u));
std_noise = zeros(1,length(bin_vec_u));
for i = 1:length(bin_vec_u)
    bin = bin_vec_u(i);
    for k = 1:K
        avg_initiation(k,i) = mean(emission_rates(k,bin_index==bin));
        std_initiation(k,i) = std(emission_rates(k,bin_index==bin));
        avg_pi0(k,i) = mean(pi0_mat(k,bin_index==bin));
        std_pi0(k,i) = std(pi0_mat(k,bin_index==bin));
        avg_occupancy(k,i) = mean(occupancy_ap(k,bin_index==bin));
        std_occupancy(k,i) = std(occupancy_ap(k,bin_index==bin));
    end
    avg_noise(i) = mean([glb_all([glb_all.bin]==bin).noise]);
    std_noise(i) = std([glb_all([glb_all.bin]==bin).noise]);
end
%Make Occupancy Fig
occ_fig = figure;
hold on
for k = 1:K
    errorbar(bin_vec_u, avg_occupancy(k,:),std_occupancy(k,:),'LineWidth',1,'Color','black');
end
plot(reshape(bin_vec_u,1,[])', avg_occupancy','black','LineWidth',1);
sctr_occ = scatter(bin_vec', reshape(occupancy_ap,1,[])' , 75, c_mat, 'o','filled', 'MarkerFaceAlpha', .2);
sctr_occ_avg = scatter(repelem(bin_vec_u,1,K), reshape(avg_occupancy,1,[]) , 75, repmat(cm,length(bin_vec_u),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',1);
box off
axis([(min(bin_range)-1) (max(bin_range)+1) 0 1.2*max(max(occupancy_ap))])
set(gca,'FontName','Lucida Sans Regular')
title('State Occupancy by AP Position');
ylabel('Occupancy Share');
xlabel('Relative Position');
grid on
saveas(occ_fig, [ OutPath '/occupancy.eps'], 'epsc');
saveas(occ_fig, [ OutPath '/occupancy.png'], 'png');

%% Make Initiation Rate Fig
emission_fig = figure;
hold on
for k = 1:K
    errorbar(bin_vec_u, avg_initiation(k,:),std_initiation(k,:),'LineWidth',1,'Color','black');
end
sctr_emission = scatter(bin_vec', init_vec , 75, c_mat, 'o','filled', 'MarkerFaceAlpha',.2) ;
sctr_emission_avg = scatter(repelem(bin_vec_u,1,K), reshape(avg_initiation,1,[]) , 75, repmat(cm,length(bin_vec_u),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',.3) ;
axis([(min(bin_range)-1) (max(bin_range)+1) 0 1.2*max(init_vec)])
set(gca,'FontName','Lucida Sans Regular')
grid on
title('State Fluorescence by AP Position');
xlabel('Relative Position');
ylabel('Fluorescence Loading Rate (A.U s^{-1})');
saveas(emission_fig, [OutPath '/emission.eps'], 'epsc');
saveas(emission_fig, [OutPath '/emission.png'], 'png');
hold off

%% Make Transition Rate Figures
%Make 3D Arrays to stack matrices
R_orig_array = zeros(K,K,length(bin_index));
R_fit_array = zeros(K,K,length(bin_index));
A_array = zeros(K,K,length(bin_index));
for i = 1:length(glb_all)
    R_fit_array(:,:,i) = reshape(glb_all(i).R_fit,K,K);
    R_orig_array(:,:,i) = real(glb_all(i).R_mat);
    A_array(:,:,i) = real(glb_all(i).A_mat);
end

% 2D Arrays to store moments (A's are calculated to have for output structure)
avg_R_orig = zeros(K,K,length(bin_vec_u));
std_R_orig = zeros(K,K,length(bin_vec_u));
avg_A = zeros(K,K,length(bin_vec_u));
std_A = zeros(K,K,length(bin_vec_u));
avg_R_fit = zeros(K,K,length(bin_vec_u));
std_R_fit = zeros(K,K,length(bin_vec_u));
for b = 1:length(bin_vec_u)
    bin = bin_vec_u(b);
    A_mean = mean(A_array(:,:,bin_index==bin),3);
    %Normalize A
    A_mean = A_mean ./ repmat(sum(A_mean),K,1);
    avg_A(:,:,b) = A_mean;    
    %Calculate and Store R moments
    avg_R_orig(:,:,b) = mean(R_orig_array(:,:,bin_index==bin),3);
    std_R_orig(:,:,b) = std(R_orig_array(:,:,bin_index==bin),1,3);
    avg_R_fit(:,:,b) = mean(R_fit_array(:,:,bin_index==bin),3);
    std_R_fit(:,:,b) = std(R_fit_array(:,:,bin_index==bin),1,3);
end

%Make multi-panel fig with fits and originals
rate_fig = figure('Position',[0 0 1024 512]);
increment = floor(60/(2*K));
for j = 1:K
    subplot(1,K,j);
    hold on
    grid on
    index = 1:K;
    index = index(index~=j);            
    legend_cell = {};
    for k = 1:(K-1)
        if plot_orig_rates            
            scatter(bin_vec_u, reshape(avg_R_orig(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'o', 'filled', 'MarkerEdgeColor', 'black') ;
        end        
        legend_cell = {legend_cell{:} ['To ' num2str(index(k))]};        
        scatter(bin_vec_u, reshape(avg_R_fit(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 's', 'filled', 'MarkerEdgeColor', 'black') ;
    end    
    legend(legend_cell{:})
    for k = 1:(K-1)
        if plot_orig_rates
            h = scatter(bin_index, reshape(R_orig_array(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'o', 'filled', 'MarkerFaceAlpha', 0.2);            
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        h = scatter(bin_index, reshape(R_fit_array(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 's', 'filled', 'MarkerFaceAlpha', 0.2) ;                
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    xlabel('Relative Position')
    ylabel('s^{-1}')    
    for k = 1:(K-1)
        if plot_orig_rates
            r_orig = reshape(avg_R_orig(index(k),j,:),1,[]);
            p = plot(bin_vec_u, r_orig,'black','LineWidth',1);
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        r_fit = reshape(avg_R_fit(index(k),j,:),1,[]);        
        p = plot(bin_vec_u, r_fit,'black','LineWidth',1);        
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    axis([min(bin_vec)-1 max(bin_vec)+1 0 1.05*max(max(max(R_fit_array)))]);
    title(['Outflow Rates from State ' num2str(j)]); 
end
saveas(rate_fig, [OutPath '/transition_rates.eps'], 'epsc');
saveas(rate_fig, [OutPath '/transition_rates.png'], 'png');
%%
rate_fig_single = figure('Position',[0 0 512 512]);
hold on
%Make Single Panel with just fits
for j = 1:K        
    grid on
    index = 1:K;
    index = index(index~=j);            
    legend_cell = {};
    for k = 1:(K-1)        
        errorbar(bin_vec_u,reshape(avg_R_fit(index(k),j,:),1,[]),reshape(std_R_fit(index(k),j,:),1,[]),'LineWidth',1,'Color','black');
        scatter(bin_index, reshape(R_fit_array(index(k),j,:),1,[]), 75, cmap(increment*((K-1)*(j-1)+k),:), 's', 'filled', 'MarkerFaceAlpha', 0.2) ;                
        scatter(bin_vec_u, reshape(avg_R_fit(index(k),j,:),1,[]), 75, cmap(increment*((K-1)*(j-1)+k),:), 's', 'filled', 'MarkerEdgeColor', 'black') ;
    end    
end
xlabel('Relative Position')
ylabel('s^{-1}')
xlim([min(bin_range) - 1 max(bin_range) + 1]) 
title('Transition Rates'); 
saveas(rate_fig_single, [OutPath '/transition_rates.eps'], 'epsc');
saveas(rate_fig_single, [OutPath '/transition_rates.png'], 'png');

%%
cycle_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(45,:)]);

bin_vec = [glb_all.bin];

off_times = dwell_times(1,:);
off_occ = occupancy_ap(1,:);
t_cycle = off_times./off_occ;
%!!! Need to update to incorporate real data
hold on

plot([glb_all.bin], t_cycle.^-1,'black','LineWidth',1);

sctr = scatter(bin_vec, t_cycle.^-1, 75, cm, 'o', 'filled', 'MarkerEdgeColor', 'black') ;
% legend('k_{off}','k_{on}');
title('Characteristic Cycle Frequency by AP Position');
xlabel('Relative Position (%)');
ylabel('Cycle Frequency (s^-1)');
grid on 
saveas(cycle_fig, [OutPath '/cycle_times.eps'], 'epsc');
saveas(cycle_fig, [OutPath '/cycle_times.png'], 'png');

n_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(5,:)]);

bin_vec = [glb_all.AP];

dp_count = [glb_all.N];
%!!! Need to update to incorporate real data
hold on

% plot([glb_all.AP], dp_count,'black','LineWidth',1);

a = area(bin_vec, dp_count);
a.FaceColor = cm;
a.EdgeColor = cm;
a.FaceAlpha = .5;
% legend('k_{off}','k_{on}');
title('Number of Data Points by AP Position');
xlabel('Relative Position (%)');
ylabel('N Data Points');
grid on 
saveas(n_fig, [OutPath '/n_dp.eps'], 'epsc');
saveas(n_fig, [OutPath '/n_dp.png'], 'png');
%%
for i = 1:length(traces_all)
    plot(traces_all(i).fluo)
    pause(1)
end

%% Make Output Struct With Relevant Fields
hmm_results = struct;
for i = 1:length(bin_vec_u)
    hmm_results(i).initiation_mean = avg_initiation(:,i);
    hmm_results(i).initiation_std = std_initiation(:,i);    
    hmm_results(i).occupancy_mean = avg_occupancy(:,i);
    hmm_results(i).occupancy_std = std_occupancy(:,i);    
    hmm_results(i).pi0_mean = avg_pi0(:,i);
    hmm_results(i).pi0_std = std_pi0(:,i);
    hmm_results(i).dwell_mean = avg_dwell(:,i);
    hmm_results(i).dwell_std = std_dwell(:,i);    
    hmm_results(i).A_mean = avg_A(:,:,i);
    hmm_results(i).R_orig_mean = avg_R_orig(:,:,i);
    hmm_results(i).R_orig_std = std_R_orig(:,:,i);
    hmm_results(i).R_fit_mean = avg_R_fit(:,:,i);
    hmm_results(i).R_fit_std = std_R_fit(:,:,i);
    hmm_results(i).noise_mean = avg_noise(i);
    hmm_results(i).noise_std = std_noise(i);
    hmm_results(i).binID = bin_vec_u(i);
    hmm_results(i).alpha = alpha;
    hmm_results(i).dT = dT;
end

save([OutPath '/hmm_results_' num2str(K) 'states_mem' num2str(w) '.mat'],'hmm_results')

