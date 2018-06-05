addpath('../utilities/');

%------------------------------Set System Params--------------------------%
%memory assumed for inference
w = 7;
%states used for final inference
K = 3;
%Tres
dT = 20;
%set bin range
bin_range = [0:7];

datatype = 'weka';
date_str = '2017-09-25';
project = ['eve7stripes_inf_2017_09_25'];

folder_path = ['../../out/' project '/' date_str '/inference_w' num2str(w) '/individual_results/'];
files = dir(folder_path);
%Include original (un-fit) transition rates in plots?
plot_orig_rates = 0;
filenames = {};
full_embryo = [];

for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
            ~isempty(strfind(files(i).name,['K' num2str(K)])) 
        filenames = [filenames {files(i).name}];
    end
    if ~isempty(strfind(files(i).name,'bin1_7'))
        full_embryo = [full_embryo 1];
    else
        full_embryo = [full_embryo 0];
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
% glb_all = glb_all(isempty([glb_all.bins_all]));

glb_all = glb_all(ismember([glb_all.bin],bin_range));
for i = 1:length(glb_all)
    if full_embryo(i) == 1
        glb_all(i).bin = 0;
    end
end

OutPath = ['../../fig/experimental_system' '/' project '/' date_str '/mem' ...
    num2str(w) '_states' num2str(K) '/' ];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

%load traces (saved as "interp_struct")
traces_all = load(['../../dat/' project '/inference_traces_t' num2str(dT) '_' project '.mat']);
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.stripe_id]),bin_range));

% %% Compare Best LogL Values across bootstraps
% for b = bin_range
%     figb = figure;
%     bin_set = glb_all([glb_all.bin]==b);
%     n_local = length([bin_set(1).local_runs]);
%     logL_scatters = zeros(n_local,length(bin_set));
%     for i = 1:length(bin_set)
%         logL_scatters(:,i) = [bin_set(i).local_runs.logL]/bin_set(i).N;
%     end
%     scatter(repelem(1:n_local,length(bin_set)),reshape(logL_scatters,1,[]),50,'MarkerFaceColor','blue','MarkerFaceAlpha',.3);
%     title(['Likelihood Scores: Stripe ' num2str(b)])
%     grid on
% end
%% Perform Rate Fitting As Necessary
%Define convenience Arrays and vectors
glb_bin_index = [glb_all.bin];
alpha = glb_all(1).alpha;
bin_index = [glb_all.bin];
bin_vec = repelem(bin_index,1,K);

%Perform Rate Fitting
dwell_times = zeros(K,length(glb_all));
for i = 1:length(glb_all)    
    [~, ranked_v] = sort(glb_all(i).v);
    A = reshape(glb_all(i).A,K,K);
    A = A(ranked_v, ranked_v);
    %Obtain raw R matrix
    R = prob_to_rate(A,dT);
    glb_all(i).R_mat = R;
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
med_dwell = zeros(K, length(bin_range));
std_dwell = zeros(K, length(bin_range));
for i = 1:length(bin_range)
    bin = bin_range(i);
    for k = 1:K
        med_dwell(k,i) = median(dwell_times(k,bin_index==bin));        
        std_dwell(k,i) = std(dwell_times(k,bin_index==bin));        
    end
end
%Make Dwell time figure
dwell_fig = figure;
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
hold on
for k = 1:K
    errorbar(bin_range, med_dwell(k,:),std_dwell(k,:),'LineWidth',1,'Color','black');
end
sctr_dwell = scatter(bin_vec', reshape(dwell_times,1,[])', 75, c_mat,...
    'o', 'filled', 'MarkerFaceAlpha', .2);
sctr_emission_avg = scatter(repelem(bin_range,1,K), reshape(med_dwell,...
    1,[]) , 75, repmat(cm,length(bin_range),1), 's','filled',...
    'MarkerEdgeColor', 'black','MarkerFaceAlpha',1);
% legend(sctr_emission_avg(:) ,legend_string{:})
set(gca,'FontName','Lucida Sans Regular')
axis([(min([glb_all.bin])-1) (max([glb_all.bin])+1) 0 1.1*max(max(dwell_times))])
title('Dwell Times by AP Position');
xlabel('Eve Stripe');
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
    A = glb_all(i).A_mat;
    A = A(ranked_r,ranked_r);
    [V,D] = eig(A);
    V = real(V);
    D = real(D);
    steady = V(:,diag(real(D))==max(diag(D)))./sum(V(:,diag(D)==max(diag(D))));
    if max(steady) > 1
        error('wtf')
    end
    occupancy_ap(:,i) = steady;
    glb_all(i).occupancy = steady;
end
init_vec = reshape(emission_rates,1,[]);

%Calculate averages
pi0_mat = vertcat(glb_all.pi0)';
med_initiation = zeros(K,length(bin_range));
std_initiation = zeros(K,length(bin_range));
med_pi0 = zeros(K,length(bin_range));
std_pi0 = zeros(K,length(bin_range));
med_occupancy = zeros(K,length(bin_range));
std_occupancy = zeros(K,length(bin_range));
med_noise = zeros(1,length(bin_range));
std_noise = zeros(1,length(bin_range));
for i = 1:length(bin_range)
    bin = bin_range(i);
    for k = 1:K
        med_initiation(k,i) = median(emission_rates(k,bin_index==bin));
        std_initiation(k,i) = std(emission_rates(k,bin_index==bin));
        med_pi0(k,i) = median(pi0_mat(k,bin_index==bin));
        std_pi0(k,i) = std(pi0_mat(k,bin_index==bin));
        med_occupancy(k,i) = median(occupancy_ap(k,bin_index==bin));
        std_occupancy(k,i) = std(occupancy_ap(k,bin_index==bin));
    end
    med_noise(i) = median([glb_all([glb_all.bin]==bin).noise]);
    std_noise(i) = std([glb_all([glb_all.bin]==bin).noise]);
end
%Make Occupancy Fig
occ_fig = figure;
hold on
for k = 1:K
    errorbar(bin_range, med_occupancy(k,:),std_occupancy(k,:),'LineWidth',1,'Color','black');
end
plot(reshape(bin_range,1,[])', med_occupancy','black','LineWidth',1);
sctr_occ = scatter(bin_vec', reshape(occupancy_ap,1,[])' , 75, c_mat, 'o','filled', 'MarkerFaceAlpha', .2);
sctr_occ_avg = scatter(repelem(bin_range,1,K), reshape(med_occupancy,1,[]) , 75, repmat(cm,length(bin_range),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',1);
box off
axis([(min(bin_range)-1) (max(bin_range)+1) 0 1.2*max(max(occupancy_ap))])
set(gca,'FontName','Lucida Sans Regular')
title('State Occupancy by AP Position');
ylabel('Occupancy Share');
xlabel('Eve Stripe');
grid on
saveas(occ_fig, [ OutPath '/occupancy.eps'], 'epsc');
saveas(occ_fig, [ OutPath '/occupancy.png'], 'png');

%% Make Initiation Rate Fig
emission_fig = figure;
hold on
for k = 1:K
    errorbar(bin_range, med_initiation(k,:),std_initiation(k,:),'LineWidth',1,'Color','black');
end
sctr_emission = scatter(bin_vec', init_vec , 75, c_mat, 'o','filled', 'MarkerFaceAlpha',.2) ;
sctr_emission_avg = scatter(repelem(bin_range,1,K), reshape(med_initiation,1,[]) , 75, repmat(cm,length(bin_range),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',.3) ;
axis([(min(bin_range)-1) (max(bin_range)+1) 1.2*abs(min(init_vec))*sign(min(init_vec)) 1.2*max(init_vec)])
set(gca,'FontName','Lucida Sans Regular')
grid on
title('State Fluorescence by AP Position');
xlabel('Eve Stripe');
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
med_R_orig = zeros(K,K,length(bin_range));
std_R_orig = zeros(K,K,length(bin_range));
med_A = zeros(K,K,length(bin_range));
std_A = zeros(K,K,length(bin_range));
med_R_fit = zeros(K,K,length(bin_range));
std_R_fit = zeros(K,K,length(bin_range));
for b = 1:length(bin_range)
    bin = bin_range(b);
    A_med = median(A_array(:,:,bin_index==bin),3);
    %Normalize A
    A_med = A_med ./ repmat(sum(A_med),K,1);
    med_A(:,:,b) = A_med;    
    %Calculate and Store R moments
    med_R_orig(:,:,b) = median(R_orig_array(:,:,bin_index==bin),3);
    std_R_orig(:,:,b) = std(R_orig_array(:,:,bin_index==bin),1,3);
    med_R_fit(:,:,b) = median(R_fit_array(:,:,bin_index==bin),3);
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
            scatter(bin_range, reshape(med_R_orig(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'o', 'filled', 'MarkerEdgeColor', 'black') ;
        end        
        legend_cell = {legend_cell{:} ['To ' num2str(index(k))]};        
        scatter(bin_range, reshape(med_R_fit(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 's', 'filled', 'MarkerEdgeColor', 'black') ;
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
    xlabel('Eve Stripe')
    ylabel('s^{-1}')
    rate_list = [];
    for k = 1:(K-1)
        if plot_orig_rates
            r_orig = reshape(med_R_orig(index(k),j,:),1,[]);
            p = plot(bin_range, r_orig,'black','LineWidth',1);
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        r_fit = reshape(med_R_fit(index(k),j,:),1,[]);        
        rate_list = [rate_list r_fit];
        p = plot(bin_range, r_fit,'black','LineWidth',1);        
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    axis([min(bin_vec)-1 max(bin_vec)+1 0 1.25*max(rate_list)]);
    title(['Outflow Rates from State ' num2str(j)]); 
end
saveas(rate_fig, [OutPath '/transition_rates.eps'], 'epsc');
saveas(rate_fig, [OutPath '/transition_rates.png'], 'png');
%% Make Output Struct With Relevant Fields
hmm_results = struct;
for i = 1:length(bin_range)
    hmm_results(i).initiation_mean = med_initiation(:,i);
    hmm_results(i).initiation_std = std_initiation(:,i);    
    hmm_results(i).occupancy_mean = med_occupancy(:,i);
    hmm_results(i).occupancy_std = std_occupancy(:,i);    
    hmm_results(i).pi0_mean = med_pi0(:,i);
    hmm_results(i).pi0_std = std_pi0(:,i);
    hmm_results(i).dwell_mean = med_dwell(:,i);
    hmm_results(i).dwell_std = std_dwell(:,i);    
    hmm_results(i).A_mean = med_A(:,:,i);
    hmm_results(i).R_orig_mean = med_R_orig(:,:,i);
    hmm_results(i).R_orig_std = std_R_orig(:,:,i);
    hmm_results(i).R_fit_mean = med_R_fit(:,:,i);
    hmm_results(i).R_fit_std = std_R_fit(:,:,i);
    hmm_results(i).noise_mean = med_noise(i);
    hmm_results(i).noise_std = std_noise(i);
    hmm_results(i).binID = bin_range(i);
    hmm_results(i).alpha = alpha;
    hmm_results(i).dT = dT;
end

save([OutPath '/hmm_results_mem' num2str(w) '_states' num2str(K) '.mat'],'hmm_results')
   
