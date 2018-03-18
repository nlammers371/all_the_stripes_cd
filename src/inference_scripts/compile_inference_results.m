% Compile results into summary structure. Generate summary plots
addpath('../utilities/');
clear all
close all
%------------------------------Set System Params--------------------------%
w = 7; %memory assumed for inference
K = 2; %states used for final inference
Tres = 20; %Time Resolution
alpha = 1.4; % MS2 rise time in time steps
fluo_type = 1; % type of spot integration used
clipped = 1; % if 0, traces are taken to be full length of nc14
stop_time_inf = 60;
clipped_ends = 1;
dynamic_bins = 1; % if 1, use time-resolved region classifications
t_window = 15; % size of sliding window used
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = stripe_range(i) - j/3 + 2/3;
    end
end
bin_range_vec = [0 bin_range_vec];
    
x_labels = {};
for i = 1:length(stripe_range)
    x_labels = [x_labels{:} {['stripe ' num2str(i) ' (A)']} {['stripe ' num2str(i) ' (C)']}...
        {['stripe ' num2str(i) ' (P)']}];
end

% id variables
datatype = 'weka';
inference_type = 'dp_bootstrap_results';
project = 'eve7stripes_inf_2018_02_20'; %project identifier


%Generate filenames and writepath
% truncated_inference_w7_t20_alpha14_f1_cl1_no_ends1
id_string = [ 'truncated_inference_w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins) '/' inference_type '/']; 
% DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\hmmm_data\inference_out\';
DropboxFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';

folder_path =  [DropboxFolder '/' project '/' id_string];
OutPath = ['../../dat/' project '/' id_string];
FigPath = ['../../fig/experimental_system/' project '/' id_string];
mkdir(OutPath)
mkdir(FigPath)
%---------------------------------Read in Files---------------------------%
files = dir(folder_path);
filenames = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
       ~isempty(strfind(files(i).name,['K' num2str(K)]))
        filenames = [filenames {files(i).name}];
    end
end

if isempty(filenames)
    error('No file with specified inference parameters found')
end

x_grid_plot = (2:22)/3;
temp_path = [DropboxFolder '/' project '/truncated_inference_w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins) '/'];
%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 1;
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'    
    load([folder_path filenames{f}]);    
    if output.skip_flag == 1 
        continue
    end    
    t_window = output.t_window;
    outpath = [temp_path '/t_window' num2str(round(t_window/60)) '/' inference_type '/'];
    mkdir(outpath)
    save([outpath filenames{f}])
    for fn = fieldnames(output)'
        glb_all(f_pass).(fn{1}) = output.(fn{1});
    end
    glb_all(f_pass).source = filenames{f};        
    f_pass = f_pass + 1
%     catch
%         warning(['File ' num2str(f) ' failed to load'])
%     end
end

%%
%%% ------------------------------Fig Calculations-------------------------%
%Adjust rates as needed (address negative off-diagonal values)
%Define convenience Arrays and vectors
alpha = glb_all(1).alpha;
bin_vec = [glb_all.stripe_id];
bin_range_vec = unique(bin_vec);
time_vec = [glb_all.t_inf];
time_index = unique(time_vec);
initiation_rates = zeros(K,length(glb_all)); % r
pi0_all = zeros(K,length(glb_all));    
noise_all = zeros(1,length(glb_all)); % sigma    
dwell_all = NaN(K,length(glb_all));
% Extract variables and perform rate fitting as needed
for i = 1:length(glb_all)    
    [initiation_rates(:,i), ranked_r] = sort(60*[glb_all(i).r]); 
    pi0 = glb_all(i).pi0(ranked_r);    
    noise_all(i) = sqrt(glb_all(i).noise); 
    A = reshape(glb_all(i).A,K,K);
    A = A(ranked_r, ranked_r);
    %Obtain raw R matrix
    R = prob_to_rate(A,Tres);
    glb_all(i).R_mat = R;
    Rcol = reshape(R,1,[]);
    glb_all(i).R = Rcol;
    %Check for imaginary and negative elements. If present, perform rate
    %fit   
    glb_all(i).r_fit_flag = 0;
    if ~isreal(Rcol)
        glb_all(i).r_fit_flag = 1;
        out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
        glb_all(i).R_fit = 60*out.R_out;
        r_diag = 60*diag(out.R_out);
    elseif sum(Rcol<0)>K
        glb_all(i).r_fit_flag = 2;
        R_conv = 60*(R - eye(K).*R);
        R_conv(R_conv<0) = 0;
        R_conv(eye(K)==1) = -sum(R_conv);
        glb_all(i).R_fit = R_conv;
        r_diag = diag(R_conv);
    else
        glb_all(i).R_fit = 60*glb_all(i).R_mat;
        r_diag = 60*diag(glb_all(i).R_mat); 
    end    
    dwell_all(:,i) = -1./r_diag;
end

% Save rate info
%Make 3D Arrays to stack matrices
R_orig_array = NaN(length(glb_all),K^2);
R_fit_array = NaN(length(glb_all),K^2);
A_array = NaN(length(glb_all),K^2);
for i = 1:length(glb_all)
    R_fit_array(i,:) = reshape(glb_all(i).R_fit,1,[]);
    R_orig_array(i,:) = reshape(real(glb_all(i).R_mat),1,[]);
    A_array(i,:) = reshape(real(glb_all(i).A_mat),1,[]);
end

% 2D Arrays to store moments (A's are calculated to have for output structure)
avg_R_orig = zeros(length(time_index),K^2,length(bin_range_vec));
std_R_orig = zeros(length(time_index),K^2,length(bin_range_vec));
avg_A = zeros(length(time_index),K^2,length(bin_range_vec));
std_A = zeros(length(time_index),K^2,length(bin_range_vec));
avg_R_fit = zeros(length(time_index),K^2,length(bin_range_vec));
std_R_fit = zeros(length(time_index),K^2,length(bin_range_vec));

for b = 1:length(bin_range_vec)
    bin = bin_range_vec(b);
    for t = 1:length(time_index)
        time = time_index(t);
        A_mean = mean(A_array(bin_vec==bin&time_vec==time,:),1);
        A_reshape = reshape(A_mean,K,K);
        %Normalize A
        A_mean = reshape(A_reshape ./ repmat(sum(A_reshape),K,1),1,[]);
        avg_A(t,:,b) = A_mean;    
        %Calculate and Store R moments
        avg_R_orig(t,:,b) = mean(R_orig_array(bin_vec==bin&time_vec==time,:));
        std_R_orig(t,:,b) = std(R_orig_array(bin_vec==bin&time_vec==time,:));
        avg_R_fit(t,:,b) = mean(R_fit_array(bin_vec==bin&time_vec==time,:));
        std_R_fit(t,:,b) = std(R_fit_array(bin_vec==bin&time_vec==time,:));
    end
end
%Occupancy
occupancy = zeros(K,length(glb_all));
for i = 1:length(glb_all)
    [~, ranked_r] = sort([glb_all(i).r]);
    A = glb_all(i).A_mat;
    A = A(ranked_r,ranked_r);
    [V,D] = eig(A);
    ind = diag(D)==max(diag(D));
    steady = V(:,ind)./sum(V(:,ind));
    occupancy(:,i) = steady;
    glb_all(i).occupancy = steady;
end
avg_dwell = NaN(K,length(bin_range_vec),length(time_index));
std_dwell = NaN(K,length(bin_range_vec),length(time_index));
avg_initiation = NaN(K,length(bin_range_vec),length(time_index));
std_initiation = NaN(K,length(bin_range_vec),length(time_index));
avg_occupancy = NaN(K,length(bin_range_vec),length(time_index));
std_occupancy = NaN(K,length(bin_range_vec),length(time_index));
avg_pi0 = NaN(K,length(bin_range_vec),length(time_index));
std_pi0 = NaN(K,length(bin_range_vec),length(time_index));
avg_noise = NaN(1,length(bin_range_vec),length(time_index));
std_noise = NaN(1,length(bin_range_vec),length(time_index));
for i = 1:length(bin_range_vec)
    bin = bin_range_vec(i);
    for t = 1:length(time_index)
        time = time_index(t);
        for k = 1:K
            avg_initiation(k,i,t) = mean(initiation_rates(k,bin_vec==bin&time_vec==time));
            std_initiation(k,i,t) = std(initiation_rates(k,bin_vec==bin&time_vec==time)); 
            avg_occupancy(k,i,t) = mean(occupancy(k,bin_vec==bin&time_vec==time));
            std_occupancy(k,i,t) = std(occupancy(k,bin_vec==bin&time_vec==time)); 
            avg_dwell(k,i,t) = mean(dwell_all(k,bin_vec==bin&time_vec==time));
            std_dwell(k,i,t) = std(dwell_all(k,bin_vec==bin&time_vec==time));        
            avg_pi0(k,i,t) = mean(pi0_all(k,bin_vec==bin&time_vec==time));
            std_pi0(k,i,t) = std(pi0_all(k,bin_vec==bin&time_vec==time));
        end  
        avg_noise(i,t) = mean(noise_all(bin_vec==bin&time_vec==time));
        std_noise(i,t) = std(noise_all(bin_vec==bin&time_vec==time));
    end
end

%%% Make Output Struct With Relevant Fields
hmm_results = struct;
for i = 1:length(bin_range_vec)
    for t = 1:length(time_index)
        ind = (i-1)*length(time_index) + t;
%         hmm_results(ind).N = bin_counts(i);
        hmm_results(ind).initiation_mean = avg_initiation(:,i,t);
        hmm_results(ind).initiation_std = std_initiation(:,i,t);    
        hmm_results(ind).occupancy_mean = avg_occupancy(:,i,t);
        hmm_results(ind).occupancy_std = std_occupancy(:,i,t);        
        hmm_results(ind).pi0_mean = avg_pi0(:,i,t);
        hmm_results(ind).pi0_std = std_pi0(:,i,t);
        hmm_results(ind).dwell_mean = avg_dwell(:,i,t);
        hmm_results(ind).dwell_std = std_dwell(:,i,t);    
        hmm_results(ind).A_mean = avg_A(t,:,i);
        hmm_results(ind).R_orig_mean = avg_R_orig(t,:,i);
        hmm_results(ind).R_orig_std = std_R_orig(t,:,i);
        hmm_results(ind).R_fit_mean = avg_R_fit(t,:,i);
        hmm_results(ind).R_fit_std = std_R_fit(t,:,i);
        hmm_results(ind).noise_mean = avg_noise(i,t);
        hmm_results(ind).noise_std = std_noise(i,t);
        hmm_results(ind).binID = bin_range_vec(i);
        hmm_results(ind).t_inf = time_index(t);
        hmm_results(ind).alpha = alpha;
        hmm_results(ind).dT = Tres;
        hmm_results(ind).fluo_type = fluo_type; 
        hmm_results(ind).clipped = clipped; 
        hmm_results(ind).clipped = clipped_ends; 
    end
end
save([OutPath '/hmm_results_mem' num2str(w)  '_states' num2str(K)  '.mat'],'hmm_results')
%%
%%% Make HMM Result Summary Figures
plot_times = time_index([2,4,6,8]);
close all
% make sub-directory for HMM plots
hmmPath = [FigPath '/hmm_figs/'];
mkdir(hmmPath);
cm = jet(128);
increment = floor(size(cm,1)/K);
sub_inc = floor(increment/length(plot_times));
color_pallette = zeros(K,3,length(plot_times));
for k = 1:K
    for t = 1:length(plot_times)
        color_pallette(k,:,t) = cm((k-1)*increment + (t-1)*sub_inc + 1,:);
    end
end


for j = 1:K
    rate_fig = figure;
    hold on    
    index = 1:K;
    index = index(index~=j);                
    for k = 1:(K-1)
        for t = 1:length(plot_times)
            col = sub2ind([K,K],index(k),j);    
            e = errorbar(bin_range_vec,reshape(avg_R_fit(time_index==plot_times(t),col,:),1,[]),...
                reshape(std_R_fit(time_index==plot_times(t),col,:),1,[]),'LineWidth',1,'Color','black');        
            e.CapSize = 0;
            if plot_scatters
                scatter(bin_vec(time_vec==plot_times(t)), R_fit_array(time_vec==plot_times(t),col), ...
                    MarkerSize, color_pallette(index(k),:,t), 's', 'filled', 'MarkerFaceAlpha', 0.2,...
                    'MarkerEdgeAlpha',0) ;                            
            end
        end
    end    
    for k = 1:(K-1)        
        for t = 1:length(plot_times)
            r_fit = reshape(avg_R_fit(time_index==plot_times(t),col,:),1,[]);        
            plot(bin_range_vec, r_fit,'black','LineWidth',1);  
        end
    end
    legend_cell = {};          
    s = [];
    for k = 1:(K-1)        
        for t = 1:length(plot_times)
            legend_cell = {legend_cell{:} ['To ' num2str(index(k)) ' (t=' num2str(round(plot_times(t)/60)) ')']};        
            s = [s scatter(bin_range_vec, reshape(avg_R_fit(time_index==plot_times(t),col,:),1,[]),...
                MarkerSize, color_pallette(index(k),:,t), 'o', 'filled', 'MarkerEdgeColor', 'black')];
        end
    end  
    legend(s,legend_cell{:},'Location','southeast')
%     axis([min(bin_range)-1 max(bin_range) + 1 0 max(max(max(R_fit_array(index,j,:))))*1.05]);
    title(['Outflow Rates from State ' num2str(j)]); 
    xlabel('stripe')
    ylabel('events per minute') 
    axis([0 8 0 max(max(avg_R_fit(:,col,:)))])
    set(gca,'xtick',1:7,'xticklabel',1:7)
    saveas(rate_fig, [hmmPath '/rates_from' num2str(j) '.png'], 'png');
    saveas(rate_fig, [hmmPath '/rates_from' num2str(j) '.pdf'], 'pdf');
end
%%
% Initiation Rates
init_fig  = figure;
hold on
% if plot_scatters
%     for k = 1:K
%         scatter(bin_vec, initiation_rates(k,:),MarkerSize,color_cell{k},'s',...
%             'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
%     end
% end
for k = 1:K
    for t = 1:length(plot_times)
        e = errorbar(bin_range_vec, avg_initiation(k,:,time_index==plot_times(t)),std_initiation(k,:,time_index==plot_times(t))...
            ,'LineWidth',1,'Color','black');
        e.CapSize = 0;
        scatter(bin_range_vec, avg_initiation(k,:,time_index==plot_times(t)), MarkerSize, color_pallette(k,:,t),...
            'o','filled', 'MarkerEdgeColor', 'black') ;
    end
end
axis([(min(bin_range_vec)-1) (max(bin_range_vec)+1) 0 1.1*max(initiation_rates(:))])
title('Initiation Rate');
xlabel('stripe');
ylabel('Initiation Rate (A.U per minute)');
set(gca,'xtick',1:7,'xticklabel',1:7)
saveas(init_fig, [hmmPath '/initiation_rates.png'], 'png');
saveas(init_fig, [hmmPath '/initition_rates.pdf'], 'pdf');
%%
% Occupancy
occ_fig = figure;
hold on

for k = 1:K
    for t = 1:length(plot_times)
        e = errorbar(bin_range_vec, avg_occupancy(k,:,time_index==plot_times(t)),...
            std_occupancy(k,:,time_index==plot_times(t)),'LineWidth',1,'Color','black');
        e.CapSize = 0;
        scatter(bin_range_vec, avg_occupancy(k,:,time_index==plot_times(t)), MarkerSize, color_pallette(k,:,t),...
            'o','filled', 'MarkerEdgeColor', 'black') ;
    end
end
axis([(min(bin_range_vec)-1) (max(bin_range_vec)+1) 0 1.2*max(max(occupancy))])
title('State Occupancy by AP Position');
ylabel('Occupancy Share');
xlabel('stripe');
set(gca,'xtick',1:7,'xticklabel',1:7)
saveas(occ_fig, [hmmPath '/occupancy_trends.png'], 'png');
saveas(occ_fig, [hmmPath '/occupancy_trends.pdf'], 'pdf');
%%
% Dwell times
dwell_fig = figure;
hold on
if plot_scatters
    for k = 1:K
        scatter(bin_vec, dwell_all(k,:),MarkerSize,color_cell{k},'s',...
            'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);
    end
end
for k = 1:K
    e = errorbar(bin_range_vec, avg_dwell(k,:),std_dwell(k,:),'LineWidth',1,'Color','black');
    e.CapSize = 0;
    scatter(bin_range_vec, avg_dwell(k,:), MarkerSize, color_cell{k},...
        'o','filled', 'MarkerEdgeColor', 'black') ;
end
axis([(min(bin_range_vec)-1) (max(bin_range_vec)+1) 0 1.2*max(avg_dwell(:))])
title('Dwell Times by AP Position');
ylabel('Dwell Times (min)');
xlabel('relative AP position (%)');
set(gca,'xtick',bin_range_vec,'xticklabel',x_labels)
saveas(dwell_fig, [hmmPath '/dwell_time_trends.png'], 'png');
saveas(dwell_fig, [hmmPath '/dwell_time_trends.pdf'], 'pdf');

% Noise
noise_fig = figure;
hold on
if plot_scatters    
    scatter(bin_vec, noise_all,MarkerSize,[.3 .3 .3],'s',...
            'filled', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0);    
end
e = errorbar(bin_range_vec, avg_noise,std_noise,'LineWidth',1,'Color',[.3 .3 .3]);
e.CapSize = 0;
scatter(bin_range_vec, avg_noise, MarkerSize, [.3 .3 .3],...
    'o','filled', 'MarkerEdgeColor', 'black') ;

axis([(min(bin_range_vec)-1) (max(bin_range_vec)+1) 0 1.2*max(avg_noise(:))])
title('Estimated Noise by AP Position');
ylabel('Noise (AU)');
xlabel('relative AP position (%)');
set(gca,'xtick',bin_range_vec,'xticklabel',x_labels)
saveas(noise_fig, [hmmPath '/noise_trends.png'], 'png');
saveas(noise_fig, [hmmPath '/noise_trends.pdf'], 'pdf');
