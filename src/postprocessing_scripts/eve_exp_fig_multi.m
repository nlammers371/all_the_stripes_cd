addpath('../utilities/');

%------------------------------Set System Params--------------------------%
%memory assumed for inference
w = 8;
%states used for final inference
K = 3;
%set bin range
bin_range = [1:2];
start_time = 0;
stop_time = 35;
datatype = 'weka';
date_str = '2017-09-08_test';
project = ['eve7stripes_inf'];
interp_data_name = ['_w' num2str(w)];
raw_data_name = [''];
folder_path = ['../../out/' project '/' date_str '/inference' interp_data_name '/individual_results/'];
files = dir(folder_path);

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
    num2str(w) '_states' num2str(K) '_s' num2str(start_time) '_' num2str(stop_time) '/' ];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

%load traces (saved as "interp_struct")
traces_all = load(['../../dat/' project '/inference_traces' interp_data_name '_' project '.mat']);
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.stripe_id]),bin_range));
% load raw traces
traces_raw = load(['../../dat/' project '/raw_traces' raw_data_name '_' project '.mat']);
traces_raw = traces_raw.trace_struct;
traces_raw = traces_raw(ismember(floor([traces_raw.binID]),bin_range));
% 
% % load ellipse info
% nuclei_all = load(['../../dat/' project '/ellipse_info' '_s' num2str(start_time) '_' num2str(stop_time) '_' project '.mat']);
% nuclei_all = nuclei_all.new_schnitz_struct;
% for i = 1:length(nuclei_all)
%     nuclei_all(i).binID = nuclei_all(i).stripe_region/2 - 3;
% end
% nuclei_all = nuclei_all(ismember(floor([nuclei_all.binID]),bin_range));
% %Remove traces from raw set for wihich we have no matching nucleus (need to
% %look into this)
% rm_set = [];
% set_vec_trace = [traces_raw.setID];
% nuc_vec_trace = [traces_raw.Nucleus];
% set_vec_nuc = [nuclei_all.setID];
% nuc_vec_nuc = [nuclei_all.Nucleus];
% for i = 1:length(traces_raw)
%     trace_nc = nuc_vec_trace(i);
%     trace_set = set_vec_trace(i);
%     
%     sets = set_vec_nuc(nuc_vec_nuc==trace_nc);
%     matches = set_vec_nuc(sets==trace_set);
%     if isempty(matches)
%         rm_set = [rm_set i];
%     elseif length(matches) > 1
%         error('degenerate identifiers');
%     end
% end
% ind_vec = 1:length(traces_raw);
% traces_raw = traces_raw(~ismember(ind_vec,rm_set));    

%% Generate Occupancy-Loading Rate Plot
%filter for desired AP range
glb_bin_index = [glb_all.bin];
alpha = glb_all(1).alpha;
occupancy_ap = zeros(K,length(glb_all));
occupancy_ap_v = zeros(K,length(glb_all));
emission_rates = zeros(K,length(glb_all));
for i = 1:length(glb_all)
    [emission_rates(:,i), ranked_r] = sort([glb_all(i).r]);
    A_log = reshape(glb_all(i).A_log,K,K);
    A_log = A_log(ranked_r, ranked_r);
    noise = glb_all(i).noise;
    pi0_log = glb_all(i).pi0_log;
    pi0_log = pi0_log(ranked_r);
                                   
    [V,D] = eig(exp(A_log));
    steady = V(:,1)./sum(V(:,1));
    occupancy_ap(:,i) = steady;
    glb_all(i).occupancy = steady;
end

%create bubble chart conveying dwell times and emission rates by AP
emission_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
if K == 3
    cm = flipud([cmap(50,:) ; cmap(30,:); cmap(20,:)]);
elseif K == 2
    cm = flipud([cmap(30,:); cmap(20,:)]);
end

init_vec = reshape(emission_rates,1,[]);
bin_vec = repelem([glb_all.bin],1,K);
c_mat = repmat(cm,length([glb_all.bin]),1);
%Calculate averages
bin_vec_u = unique([glb_all.bin]);
bin_index = [glb_all.bin];
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

%calculate averages
hold on

plot(reshape(bin_vec_u,1,[])', avg_initiation','black','LineWidth',1);
sctr_emission = scatter(bin_vec', init_vec , 75, c_mat, 'o','filled', 'MarkerFaceAlpha',.2) ;
sctr_emission_avg = scatter(repelem(bin_vec_u,1,K), reshape(avg_initiation,1,[]) , 75, repmat(cm,length(bin_vec_u),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',.3) ;
axis([(min(bin_range)-1) (max(bin_range)+1) 0 1.2*max(init_vec)])

set(gca,'FontName','Lucida Sans Regular')
grid on
title('State Fluorescence by AP Position');
xlabel('Relative Position');
ylabel('Fluorescence Loading Rate (A.U s^{-1})');
% ylabel(h, 'Average Occupancy')

saveas(emission_fig, [OutPath '/emission.eps'], 'epsc');
saveas(emission_fig, [OutPath '/emission.png'], 'png');
hold off

%Occupancy
occ_fig = figure;
hold on
plot(reshape(bin_vec_u,1,[])', avg_occupancy','black','LineWidth',1);
sctr_occ = scatter(bin_vec', reshape(occupancy_ap,1,[])' , 75, c_mat, 'o','filled', 'MarkerFaceAlpha', .2);
sctr_occ_avg = scatter(repelem(bin_vec_u,1,K), reshape(avg_occupancy,1,[]) , 75, repmat(cm,length(bin_vec_u),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',1);
% h;
box off
axis([(min(bin_range)-1) (max(bin_range)+1) 0 1.2*max(max(occupancy_ap))])

set(gca,'FontName','Lucida Sans Regular')

title('State Occupancy by AP Position');
ylabel('Occupancy Share');
xlabel('Relative Position');

grid on

saveas(occ_fig, [ OutPath '/occupancy.eps'], 'epsc');
saveas(occ_fig, [ OutPath '/occupancy.png'], 'png');

%% Dwell Times
dwell_times = zeros(K,length(glb_all));
dT = 22.5;%traces_all(1).dT;
for i = 1:length(glb_all)
    Rcol = glb_all(i).R;
%     R = reshape(Rcol,K,K);
    [~, ranked_v] = sort(glb_all(i).v);
    A = reshape(glb_all(i).A,K,K);
    A = A(ranked_v, ranked_v);
    R = reshape(glb_all(i).R,K,K);
    glb_all(i).R_mat = R(ranked_v, ranked_v);
    glb_all(i).R = reshape(R(ranked_v, ranked_v),1,[]);
   
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
%%
%Calculate average dwell times
avg_dwell = zeros(K, length(bin_vec_u));
std_dwell = zeros(K, length(bin_vec_u));
for i = 1:length(bin_vec_u)
    bin = bin_vec_u(i);
    for k = 1:K
        avg_dwell(k,i) = mean(dwell_times(k,bin_index==bin));        
        std_dwell(k,i) = std(dwell_times(k,bin_index==bin));        
    end
end

dwell_fig = figure;
hold on
plot(bin_vec_u', avg_dwell','black','LineWidth',1);

sctr_dwell = scatter(bin_vec', reshape(dwell_times,1,[])', 75, c_mat, 'o', 'filled', 'MarkerFaceAlpha', .2);
sctr_emission_avg = scatter(repelem(bin_vec_u,1,K), reshape(avg_dwell,1,[]) , 75, repmat(cm,length(bin_vec_u),1), 's','filled', 'MarkerEdgeColor', 'black','MarkerFaceAlpha',.3) ;
% legend(h, 'OFF','ON1','ON2');

box off
set(gca,'FontName','Lucida Sans Regular')
axis([(min([glb_all.bin])-1) (max([glb_all.bin])+1) 0 15*ceil(max(max(dwell_times))/15)])
title('Dwell Times by AP Position');
xlabel('Relative Position');
ylabel('Dwell Time (s)');
grid on 
saveas(dwell_fig, [OutPath '/dwell_times.eps'], 'epsc');
saveas(dwell_fig, [OutPath '/dwell_times.png'], 'png');

%% Look at Rates
rate_fig = figure('Position',[0 0 1024 512]);

colormap((colormap('jet')));
cmap = colormap;
bin_vec = [glb_all.bin];
R_orig_array = zeros(K,K,length(bin_vec));
R_fit_array = zeros(K,K,length(bin_vec));
A_array = zeros(K,K,length(bin_vec));
for i = 1:length(glb_all)
    R_fit_array(:,:,i) = reshape(glb_all(i).R_fit,K,K);
    R_orig_array(:,:,i) = real(glb_all(i).R_mat);
    A_array(:,:,i) = real(glb_all(i).A_mat);
end
avg_R_orig = zeros(K,K,length(bin_vec_u));
std_R_orig = zeros(K,K,length(bin_vec_u));
avg_A = zeros(K,K,length(bin_vec_u));
std_A = zeros(K,K,length(bin_vec_u));
avg_R_fit = zeros(K,K,length(bin_vec_u));
std_R_fit = zeros(K,K,length(bin_vec_u));
for b = 1:length(bin_vec_u)
    bin = bin_vec_u(b);
    A_mean = mean(A_array(:,:,bin_vec==bin),3);
    A_mean = A_mean ./ repmat(sum(A_mean),K,1);
    avg_A(:,:,b) = A_mean;    
    avg_R_orig(:,:,b) = mean(R_orig_array(:,:,bin_vec==bin),3);
    std_R_orig(:,:,b) = std(R_orig_array(:,:,bin_vec==bin),1,3);
    avg_R_fit(:,:,b) = mean(R_fit_array(:,:,bin_vec==bin),3);
    std_R_fit(:,:,b) = std(R_fit_array(:,:,bin_vec==bin),1,3);
end

increment = floor(60/(2*K));
for j = 1:K
    subplot(1,K,j);
    hold on
    grid on
    index = 1:K;
    index = index(index~=j);            
    legend_cell = {};
    for k = 1:(K-1)
        scatter(bin_vec, reshape(R_orig_array(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'o', 'filled', 'MarkerFaceAlpha', 0.2) ;
        scatter(bin_vec, reshape(R_fit_array(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'd', 'filled', 'MarkerFaceAlpha', 0.2) ;
        legend_cell = {legend_cell{:} ['R' num2str(index(k)) num2str(j) ' Original'] ['R' num2str(index(k)) num2str(j) ' Fit']};
        scatter(bin_vec_u, reshape(avg_R_orig(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'o', 'filled', 'MarkerEdgeColor', 'black') ;
        scatter(bin_vec_u, reshape(avg_R_fit(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 'd', 'filled', 'MarkerEdgeColor', 'black') ;
    end
%     axis([(min(ap_vec)-1) (max(ap_vec)+1) -.03 1.2*max(max(max(([R_orig_array R_fit_array]))))])
    xlabel('Relative Position')
    ylabel('s^{-1}')
%     legend(legend_cell{:},'Location', 'northeast');
    ymax = -Inf;
    ymin = 0;
    for k = 1:(K-1)
        r_orig = reshape(avg_R_orig(index(k),j,:),1,[]);
        r_fit = reshape(avg_R_fit(index(k),j,:),1,[]);
        ymax = max([r_orig,r_fit,ymax]);
        ymin = min([r_orig,r_fit,ymin]);
        plot(bin_vec_u, r_orig,'black','LineWidth',1);
        plot(bin_vec_u, r_fit,'black','LineWidth',1);        
    end
    axis([min(bin_vec)-1 max(bin_vec) + 1 1.1*ymin 1.3*ymax]);
    title(['Outflow Rates from State ' num2str(j)]); 
end
saveas(rate_fig, [OutPath '/transition_rates.eps'], 'epsc');
saveas(rate_fig, [OutPath '/transition_rates.png'], 'png');

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

%% Mean and Cumulative Trace Profiles

%Plot Mean Fluorescence as a Function of AP Region
mean_fluo_time = zeros(1,length(bin_range));
mean_fluo_time_clean = zeros(1,length(bin_range));
mean_fluo_nucleus = zeros(1,length(bin_range));
total_fluo = zeros(1,length(bin_range));
mean_duration = zeros(1,length(bin_range));
total_duration = zeros(1,length(bin_range));
duration_vec = [];
bin_vec = [];
for b = 1:length(bin_range)
    bin_traces_clean = traces_all([traces_all.binID]==bin_range(b));
    bin_traces = traces_raw([traces_raw.binID]==bin_range(b));
    mean_fluo_time(b) = mean([bin_traces.fluo]);
    mean_fluo_time_clean(b) = mean([bin_traces_clean.fluo]);
    f_sum = sum([bin_traces.fluo]);
    N = length(bin_traces);
    total_fluo(b) = f_sum;
    dur = zeros(1,N);
    for d = 1:N
        dur(d) = length(bin_traces(d).fluo);
%         dur(d) = max(bin_traces(d).time_full) - min(bin_traces(d).time_full);
    end
    duration_vec = [duration_vec dur];
    bin_vec = [bin_vec ones(1,length(dur))*bin_range(b)];
    mean_duration(b) = mean(dur);
    total_duration(b) = sum(dur);
    mean_fluo_nucleus(b) = f_sum / N;
end
mean_expected_productivity = sum(avg_occupancy .* avg_initiation)*(w-alpha/2)*traces_all(1).dT;
expected_productivity = sum(emission_rates.*occupancy_ap)*(w-alpha/2)*traces_all(1).dT;
on_share = zeros(1,length(bin_range));
on_count = zeros(1,length(bin_range));
for b = 1:length(bin_range)
    bin = bin_range(b);
    bin_nuclei = nuclei_all([nuclei_all.binID] == bin);
    on_share(b) = length(bin_nuclei([bin_nuclei.cm_fluo]>0)) / length(bin_nuclei);
    on_count(b) = length(bin_nuclei([bin_nuclei.cm_fluo]>0));
end
colormap('jet')
cm = colormap;
%Plot Time Average of Productivity vs. HMM Expectation
model_vs_exp = figure;
hold on
plot(bin_range,mean_expected_productivity,'-','LineWidth',2,'Color', cm(5,:))
scatter(glb_bin_index, expected_productivity,40,'MarkerFaceColor',cm(5,:),...
    'MarkerFaceAlpha',.3,'MarkerEdgeColor',cm(5,:));
plot(bin_range,mean_fluo_time_clean,'-','LineWidth',2,'Color', cm(30,:))
grid on
legend('HMM Prediction (Mean)', 'Individual HMM Results', 'Experimental Data', ...
    'Location','southwest');
title('Predicted vs Actual Average Fluorescence')
saveas(model_vs_exp, [OutPath '/mean_prod_fig.png'],'png')

making_a_stripe = figure('Position', [0 0 1024 1024]);
subplot(2,2,1);
hold on
area(bin_range,total_fluo/min(total_fluo),'FaceColor','black','FaceAlpha',.3)
area(bin_range,on_share/on_share(end),'FaceColor',cm(1,:),'FaceAlpha',.5)
legend('Full Profile', 'Fraction On Effect') 
grid on;
title(['Impact of Fraction Nuclei On Through ' num2str(stop_time) ' Minutes'])

subplot(2,2,2)
hold on
area(bin_range,total_fluo/min(total_fluo),'FaceColor','black','FaceAlpha',.3)
area(bin_range,(mean_duration)/(mean_duration(end)),'FaceColor',cm(15,:),'FaceAlpha',.5)
legend('Full Profile','Duration Effect');
grid on
title(['Impact of Activity Duration Through ' num2str(stop_time) ' Minutes'])

subplot(2,2,3);
hold on
area(bin_range,total_fluo/min(total_fluo),'FaceColor','black','FaceAlpha',.3)
area(bin_range,mean_fluo_time/mean_fluo_time(end),'FaceColor',cm(30,:),'FaceAlpha',.5)
legend('Full Profile', 'Dynamics Effect') 
grid on;
title(['Impact of Dynamics Through ' num2str(stop_time) ' Minutes'])

subplot(2,2,4);
hold on
% area(bin_range,total_fluo/min(total_fluo),'FaceColor','black','FaceAlpha',.3)
area(bin_range,(on_count.*mean_duration.*mean_fluo_time)/(on_count(end)*mean_duration(end).*mean_fluo_time(end)),'FaceColor','black','FaceAlpha',.3)
area(bin_range,(on_count.*mean_duration.*mean_fluo_time)/(on_count(end)*mean_duration(end).*mean_fluo_time(end)),'FaceColor',cm(30,:),'FaceAlpha',.5)
area(bin_range,(on_count.*mean_duration)/(on_count(end)*mean_duration(end)),'FaceColor',cm(15,:),'FaceAlpha',.5)
area(bin_range,on_share/on_share(end),'FaceColor',cm(1,:),'FaceAlpha',.3)
grid on
% legend('Fraction On', 'Fraction On + Duration', 'Fraction On + Duration + Differential Productivity')
title(['Combined Effects Through ' num2str(stop_time) ' Minutes'])
saveas(making_a_stripe, [OutPath '/stripe_factors_fig.png'],'png')

%Scatter plot of trace lengths
duration_fig = figure;
hold on
scatter(bin_vec,duration_vec,30,'MarkerFaceColor',cm(20,:),'MarkerFaceAlpha',.05,'MarkerEdgeColor',cm(20,:),'MarkerEdgeAlpha',.05);
scatter(bin_range,mean_duration,40,'MarkerFaceColor',cm(20,:),'MarkerFaceAlpha',1,'MarkerEdgeColor','black');
axis([min(bin_range)-1 max(bin_range) + 1 0 max(duration_vec)*1.05])
title('Trace Duration by Region')
grid on
saveas(duration_fig, [OutPath '/trace_lengths_fig.png'],'png')

%Make plots of activity by nucleus
nuclei_prod_fig = figure;
set_list = unique([nuclei_all.setID]);
nc_plot_path = [OutPath '/nc_plots/'];
if exist(nc_plot_path) ~= 7
    mkdir(nc_plot_path)
end
for s = 1:length(set_list)
    nc_fig = figure;
    set_num = set_list(s);
    set_nuclei = nuclei_all([nuclei_all.setID]==set_num);
    x_vec = [set_nuclei.xMean];
    y_vec = [set_nuclei.yMean];
    f_vec = [set_nuclei.cm_fluo];
    sc = scatter(x_vec, y_vec, 70, f_vec, 'filled');
    set(sc,'MarkerFaceAlpha',.5)
    grid on
    title(['Cumulative Productivity by Nucleus Set ' num2str(set_num)...
        ' (Mean Positions Used)'])
    saveas(nc_fig, [nc_plot_path '/nc_prod_fig_set_' ...
        num2str(set_num) '.png'],'png')
end