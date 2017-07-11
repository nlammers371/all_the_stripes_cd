%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'C:\Users\Nicholas\Dropbox (Garcia Lab)\DropboxSingleTraces\Eve2_ML\';
% folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\DropboxSingleTraces\Eve2_ML';
project = 'mHMMeve2_ml_testing';
% outName = 'eve2Sets_2017_06_15_ml.mat'
outpath = ['../projects/' project '/' ];
% Keyword to ensure only sets from current project are pulled
keyword = 'eve2_20sec_';
exclude = 'eve2_20sec_5';


%-------------------------Set Summary Parameters--------------------------%
nuclear_cycles = [14];
% Minimumt number of data points per summary stat
min_stat = 500;


if exist(outpath) ~= 7
    mkdir(outpath);
end
%NL: this could be made much better...
dir_struct = struct;
i_pass = 1;
dirinfo = dir(folder_path);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    % skip files lacking project keyword or containing names of skip sets
    if isempty(strfind(thisdir,keyword)) || ~isempty(strfind(thisdir,exclude))
        continue
    end
    subdir_struct = dir(fullfile(folder_path,thisdir, '*.mat'));
    dir_struct(i_pass).files = subdir_struct;
    dir_struct(i_pass).folder = thisdir;
    i_pass = i_pass + 1;
end

filenames = {};
for i = 1:length(dir_struct)
    subdir_struct = dir_struct(i).files;
    for j = 1:length(subdir_struct)
        if strfind(subdir_struct(j).name,'CompiledParticles') > 0
            filenames = [filenames {[folder_path dir_struct(i).folder '\' subdir_struct(j).name]}];
        end
    end
end
%%
%Data structure to store extracted trace sets
traces_interp = struct;
diff_traces = 0;
trace_struct = struct;
i_iter = 1;
for k = 1:length(filenames)
    raw_data  = load([filenames{k}]);
    time = raw_data.ElapsedTime*60;
%     Normalize to Start of 14
    
    traces = raw_data.AllTracesVector;
    filter = sum(raw_data.ncFilter(:,ismember([raw_data.ncFilterID],nuclear_cycles)),2)>0;
    traces = traces(:,filter);
    no_cycle_start = ~isnan(max(max(traces(time < 2*60,:)))) || max(max(traces(time < 2*60,:))) > 0;
    time = time - time(raw_data.nc14);
    for i = 1:size(traces,2)
        raw_trace = traces(:,i);
        if length(raw_trace(~isnan(raw_trace))) < 2
            continue
        end
        start = find(~isnan(raw_trace),1);
        stop = find(~isnan(raw_trace),1,'last');
        trunc_trace = [raw_trace(start:stop)'];
        trunc_time = time(start:stop); 
        [~, apPos] = max(raw_data.APFilter(i,:));
        trace_struct(i_iter).fluo = trunc_trace;
        trace_struct(i_iter).time = trunc_time;
        trace_struct(i_iter).AP = apPos;
        trace_struct(i_iter).set = filenames{k};
        trace_struct(i_iter).setID = k;
        trace_struct(i_iter).background = raw_data.MeanOffsetVector;
        trace_struct(i_iter).no_nc_start = no_cycle_start;
        i_iter = i_iter + 1;
    end
end
%% Interpolate Traces
% Run data quality checks. Discard suspect dps. Will be replaced in interp step
% Get 95th percentile for point-to-point deltas
% diffs = abs(diff([trace_struct.fluo]));
% ref_len = prctile(diffs(~isnan(diffs)),85);
% adjustments = 0;
% for i = 1:length(trace_struct)
%     trace = trace_struct(i).fluo;
%     trace(isnan(trace)) = 0;
%     tr_d = [0 diff(trace)];
%     tr_dd = [0 diff(diff(trace)) 0];
%     rm_list = [];
%     for j = 1:length(trace)
%         % remove "large" drops to zero
%         if abs(tr_d(j)) > ref_len && trace(j) == 0
%             rm_list = [rm_list j];
%             adjustments = adjustments +1;
%         % remove "large" transient spikes
%         elseif abs(tr_dd) > 1.5*ref_len
%             rm_list = [rm_list j];
%             adjustments = adjustments +1;
%         end
%     end
%     trace_struct(i).fluo = trace(~ismember(1:length(trace), rm_list));
%     trace_struct(i).time = trace_struct(i).time(~ismember(1:length(trace),rm_list));
% end
       
% Interpolate to Achieve Desired Memory
% Define Desired res and associated param
empirical_res = [];
for t = 1:length(trace_struct)
    empirical_res = [empirical_res diff([trace_struct(t).time])];
end
% p5 = prctile(emprical_res,5);
empirical_res = empirical_res(~isnan(empirical_res));
p95 = prctile(empirical_res,95);
empirical_res = empirical_res(empirical_res < p95);
mean_empirical_res = mean(empirical_res);    
T_elong = 160;
mem = round(T_elong/mean_empirical_res);
Tres = T_elong / mem;


%Set minimum trace length (in time steps)
min_len = 0;
i_pass = 1;
% define time grid. All traces will be shifted to conform to same 
% Assume that anything starting <5 min is not taken from a set with
% no nc14 normaliz
t_start = Inf;
for i = 1:length(trace_struct)
    if trace_struct(i).no_nc_start ~= 1
        t_start  = min(t_start,floor(min(trace_struct(i).time)));
    end
end

t_stop = floor(max([trace_struct.time]));
time_grid = t_start + Tres*(0:floor(t_stop/Tres));
interp_struct = struct;
for i = 1:length(trace_struct)
    fluo = trace_struct(i).fluo;
    time = trace_struct(i).time;
    
    time_interp = time_grid(time_grid >= min(time));
    time_interp = time_interp(time_interp <= max(time));
    fluo_interp = interp1(time,fluo,time_interp);
    if length(time_interp) > min_len
        interp_struct(i_pass).fluo = fluo_interp;
        interp_struct(i_pass).time = time_interp;
        interp_struct(i_pass).fluo_orig = fluo;
        interp_struct(i_pass).time_orig = time;
        interp_struct(i_pass).AP = trace_struct(i).AP;
        interp_struct(i_pass).set = trace_struct(i).set;
        interp_struct(i_pass).setID = trace_struct(i).setID;
        interp_struct(i_pass).N = length(fluo_interp);
        interp_struct(i_pass).dT = Tres;
        interp_struct(i_pass).T_elong = T_elong;
        interp_struct(i_pass).w = mem;
        interp_struct(i_pass).alpha = (60/204)*mem;
        i_pass = i_pass + 1;
    end
end
%%
for i = 1:length(interp_struct)
    hold off
   
    plot(interp_struct(i).time_orig, interp_struct(i).fluo_orig, 'Linewidth',1.5)
     hold on
%     t_fig = figure('Visible','off');
    plot(interp_struct(i).time, interp_struct(i).fluo,'-o', 'Linewidth',1.5)
    title(['Trace ' num2str(i)])
    grid on
    pause(1.5)
    
%     saveas(t_fig, [outpath,'/traces/', 'trace_' num2str(i) '.png'],'png');
end
%% Histogram of Data Points by AP

sets = unique({trace_struct.set});
n_sets = length(sets);

set_titles = {};
for i = 1:length(sets)
    start = strfind(sets{i},'\2017') + 12;
    stop = strfind(sets{i},'\Compiled') - 1;
    string = sets{i};
    set_titles = {set_titles{:}  strrep(string(start:stop),'_',' ')};        
end

figure(1);
subplot(2,1,1);

colormap('jet');
cm = colormap;

ap_vec = unique([trace_struct.AP]);
ap_ct_vec = zeros(1,length(ap_vec));
for a = 1:length(ap_vec)
    ap = ap_vec(a);
    ap_ct_vec(a) = sum(length([trace_struct([trace_struct.AP] == ap).fluo]));
end
bar(ap_vec, ap_ct_vec , 'Facecolor',cm(5,:),...
    'EdgeColor','black','FaceAlpha',.65,'BarWidth', 1);
title('Data Points by AP Region');
xlabel('AP Region (%)');

subplot(2,1,2);
histogram([trace_struct.AP],'Facecolor',cm(5,:));
title('Traces by AP Region');
xlabel('AP Region (%)');

saveas(figure(1), [outpath, 'data_histograms.png'],'png');

%% Run Quality Checks by data set and AP
ap_vec = unique([trace_struct.AP]);
% Determine appropriate AP grouping
if min(ap_ct_vec) >= min_stat
    ap_grp_vec = ap_vec;
else
    for i = 2:length(ap_vec)
        dupe_flag = 0;
        ap_grp_ind = floor(ap_vec / i);
        ap_grps = unique(ap_grp_ind);
        for g = ap_grps
            ct = ap_ct_vec(ap_grp_ind==g);
            if ct < min_stat
                dupe_flag = dupe_flag + 1;
                break
            end
        end
        if dupe_flag == 0
            grp = i;
            break
        end
    end
end
        
ap_names = {};
for ap_grp = ap_grps
    ap_names = {ap_names{:} [num2str(min(ap_vec(ap_grp_ind==ap_grp))) '-' num2str(max(ap_vec(ap_grp_ind==ap_grp)))]};
end  
    
n_grps = length(ap_grps);

increment = floor(60 / (n_grps *n_sets));
granularity = 100;
max_fluo = ceil(max([trace_struct.fluo])/granularity)*granularity;
FluoBins = 1:20:max_fluo;

figure(2);
for i = 1:n_grps
    ap_struct = trace_struct(ismember([trace_struct.AP],ap_vec(ap_grp_ind==ap_grps(i))));
    for j = 1:n_sets
        f_list = [];
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                fluo = ap_struct(a).fluo;
                f_list = [f_list fluo(fluo>0)];
            end
        end
        if isempty(f_list)
            continue
        end
        ap_ct = histc(f_list, FluoBins);        
        subplot(n_sets,n_grps, (j-1)*n_grps + i)
        bar(FluoBins, ap_ct / max(ap_ct), 'FaceColor',cm(increment*((j-1)*n_grps + i),:),...
            'EdgeColor',cm(increment*((j-1)*n_grps + i),:));
        set(gca,'fontsize',4)
        set_start = strfind(sets{j},'_');
        set_start = set_start(end) + 1;
        set_end = strfind(sets{j},'\Comp') - 1;
        set_num = sets{j}(set_start:set_end);
        title(['Fluo, Set: '  num2str(set_num) ' AP: ' ap_names{i}]); %' Set:' sets{j}])
        axis([0,max_fluo,0 ,1])
    end
end
hold off
saveas(figure(2), [outpath, 'fluo_his.png'],'png');
%% Cumulative Fluorescence within Eve Stripe 2 Region
%Make Strings for legen entries
ap_range = 39:45;
figure(3)
increment = floor(60 / n_sets);
ptile_list = zeros(1,n_sets);
% Set scale for x axis. Max of selected percentile across sets will be used
ptile = 97;
hold on
for j = 1:n_sets
    f_list = [];
    for a = 1:length(trace_struct)
        if strcmp(trace_struct(a).set,sets{j}) && ismember(trace_struct(a).AP,ap_range)
            f_list = [f_list trace_struct(a).fluo];
        end
    end
    ptile_list(j) = prctile(f_list,ptile); 
    if isempty(f_list)
        continue
    end
    ap_ct = histc(f_list, FluoBins);
    plot(FluoBins, cumsum(ap_ct) / sum(ap_ct), 'Color',cm(increment*j,:),'LineWidth', 2);
    
end 
title(['Cumulative PDF (AP: ' num2str(min(ap_range)) '-' num2str(max(ap_range)) ')']); 
axis([0,max(ptile_list),0,1])
grid on
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(figure(3), [outpath, 'cum_his.png'],'png');

%% Plot Fluo Percentiles By AP and Data Set

med_fig = figure('Position',[0 0 1024 512]);

medians = nan(n_sets,length(ap_vec));
nTile25 = nan(n_sets,length(ap_vec));
nTile75 = nan(n_sets,length(ap_vec));
% hold on
for j = 1:n_sets
    for k = 1:length(ap_vec)
        f_list = [];
        ap_grp = ap_vec(max(1,k-1):min(length(ap_vec),k+1));
        ap_struct = trace_struct(ismember([trace_struct.AP],ap_grp));
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                
                f_list = [f_list ap_struct(a).fluo];
            end
        end
        if isempty(f_list)
            continue
        end
        f_list = f_list(f_list > 0);
%         f_list = f_list(f_list < 1500);
        medians(j,k) = median(f_list);
        nTile25(j,k) = prctile(f_list,25);
        nTile75(j,k) = prctile(f_list,75);
    end 
end

hold on

for i = 1:n_sets
    p25_vec = nTile25(i,:); 
    p75_vec = nTile75(i,:); 
    p=[p25_vec(~isnan(p25_vec)),fliplr(p75_vec(~isnan(p25_vec)))];  
    x = [ap_vec(~isnan(p25_vec)) fliplr(ap_vec(~isnan(p25_vec)))];         
    h = fill(x,p,cm(increment*i,:),'Linewidth',1,'Edgecolor',cm(increment*i,:));  
    set(h,'facealpha',.05)
    set(h,'edgealpha',.3)
end
for i = 1:n_sets
    plot(ap_vec,medians(i,:), 'Color', cm(increment*i,:),'Linewidth',2)
    legend(set_titles{:}, 'Location','southeast')
end
title('Median Fluorescent Intensities by AP Position and Data Set')
xlabel('AP Position (%)')
ylabel('AU');
grid on
hold off
saveas(med_fig, [outpath, 'median_fluo_plots.png'],'png');

%% Mean Fluo By Region and Set
flank_anterior = 33:38;
flank_posterior = 46:51;
stripe2 = 39:45;
grp_names = {'Anterior Flank', 'Eve2 Stripe', 'Posterior Flank'};
grp_indices = {flank_anterior, stripe2, flank_posterior};
increment = floor(60 /n_sets);
tr_gross = 0;

for k = 1:length(grp_indices)
    legend_names = {};
    fig =  figure('Position',[0 0 1024 512]);
    hold on
    ap_grp = grp_indices{k};
    ap_struct = trace_struct(ismember([trace_struct.AP],ap_grp));
     
    for j = 1:n_sets
        set_traces = ap_struct([ap_struct.setID]==j);
        if isempty(set_traces)
            continue
        end
        if set_traces(1).no_nc_start == 1
            continue
        end
        legend_names = [legend_names set_titles{j}];
        set = set_traces(1).set;
        times = [set_traces.time];         
        fluo_values = [set_traces.fluo];
        fluo_values = fluo_values(~isnan(fluo_values));
        times = times(~isnan(fluo_values));
        unique_times = sort(unique(times));
        f_series = zeros(1,length(unique_times));
        for t = 1:length(unique_times)
            f_series(t) = mean(fluo_values(times == unique_times(t)));
        end
        plot(unique_times, f_series, 'Color', cm(j*increment,:),'Linewidth',2);
    end
    title(['Mean Fluorescence by Data Set: ' grp_names{k}]) %' (AP: ' num2str(ap_grp) ')'])
    xlabel('seconds');
    ylabel('AU');
    
    grid on;
    legend(legend_names{:}, 'Location','southwest')
    saveas(fig, [outpath, 'mean_temp_fluo_' grp_names{k} '.png'],'png');
    
end



