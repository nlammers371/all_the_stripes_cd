%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\DropboxSingleTraces\Eve2_ML';
project = 'mHMMeve2_ml';
% outName = 'eve2Sets_2017_06_15_ml.mat';
outpath = ['../projects/' project '/' ];
% Keyword to ensure only sets from current project are pulled
keyword = 'eve2_20sec_';
exclude = 'eve2_20sec_5';
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
  % skip files lacking prject keyword
  if isempty(strfind(thisdir,keyword)) || ~isempty(strfind(thisdir,exclude))
      continue
  end
  subdir_struct = dir(fullfile(folder_path,thisdir, '*.mat'));
  if length(subdir_struct) > 5
      if i_pass == 1
          dir_struct(i_pass).files = subdir_struct;
          i_pass = 2;
      else
          dir_struct(i_pass).files = subdir_struct;
          i_pass = i_pass + 1;
      end
  end
end

filenames = {};
for i = 1:length(dir_struct)
    subdir_struct = dir_struct(i).files;
    for j = 1:length(subdir_struct)
        if strfind(subdir_struct(j).name,'CompiledParticles') > 0
            filenames = [filenames {[subdir_struct(j).folder '\' subdir_struct(j).name]}];
        end
    end
end
%%
traces_no_nan = struct;
traces_interp = struct;
diff_traces = 0;
%Data structure to store extracted trace sets
trace_struct = struct;
i_iter = 1;
for k = 1:length(filenames)
    raw_data  = load([filenames{k}]);
    time = raw_data.ElapsedTime*60;
    time = time - time(raw_data.nc14);
    traces = raw_data.AllTracesVector;
    nc14_filter = raw_data.ncFilter(:,end);
    traces_14 = traces(:,nc14_filter == 1);
    %Set Nans to 0, true zeros to -1
    traces_14(traces_14==0) = -1000;
    traces_14(isnan(traces_14)) = 0;
    for i = 1:size(traces_14,2)
        raw_trace = traces_14(:,i);
        % Skip single points
        if length(raw_trace(raw_trace ~= 0)) < 2
            continue
        end
        start = find(raw_trace,1);
        stop = find(raw_trace,1,'last');
        raw_trace(raw_trace==0) = nan;
        raw_trace(raw_trace==-1000) = 0;
        trunc_trace = [raw_trace(start:stop)'];
        trunc_time = time(start:stop); 
        [~, apPos] = max(raw_data.APFilter(i,:));
        
        trace_struct(i_iter).fluo = trunc_trace;
        trace_struct(i_iter).time = trunc_time;
        trace_struct(i_iter).AP = apPos;
        trace_struct(i_iter).set = filenames{k};
        
        i_iter = i_iter + 1;
    end
end

% Run data quality checks. Discard suspect dps. Will be replaced in interp step
% Get 95th percentile for point-to-point deltas
diffs = abs(diff([trace_struct.fluo]));
ref_len = prctile(diffs(~isnan(diffs)),85);
adjustments = 0;
for i = 1:length(trace_struct)
    trace = trace_struct(i).fluo;
    trace(isnan(trace)) = 0;
    tr_d = [0 diff(trace)];
    tr_dd = [0 diff(diff(trace)) 0];
    rm_list = [];
    for j = 1:length(trace)
        % remove "large" drops to zero
        if abs(tr_d(j)) > ref_len && trace(j) == 0
            rm_list = [rm_list j];
            adjustments = adjustments +1;
        % remove "large" transient spikes
        elseif abs(tr_dd) > 1.5*ref_len
            rm_list = [rm_list j];
            adjustments = adjustments +1;
        end
    end
    trace_struct(i).fluo = trace(~ismember(1:length(trace), rm_list));
    trace_struct(i).time = trace_struct(i).time(~ismember(1:length(trace),rm_list));
end
       
% Interpolate to Achieve Desired Memory
% Define Desired res and associated param
mem = 20;
T_elong = 160;
Tres = T_elong / mem;
%Set minimum trace length (in time steps)
min_len = mem;

i_pass = 1;
% define time grid. All traces will be shifted to conform to same grid
t_start = Inf;
for i = 1:length(trace_struct)
    if isempty(strfind(trace_struct(i).set,'eve2_20sec_3'))
        t_start  = min(t_start,floor(min(trace_struct(i).time)));
    end
end
t_stop = floor(max([trace_struct.time]));
time_grid = t_start + Tres*(0:floor(t_stop/Tres));
interp_struct = struct;
for i = 1:length(trace_struct)
    if ~isempty(strfind(trace_struct(i).set,'eve2_20sec_3'))
        offset = t_start;
    else
        offset =0;
    end
    fluo = trace_struct(i).fluo;
    time = offset + trace_struct(i).time;
    
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
%     plot(interp_struct(i).time, interp_struct(i).fluo, 'Linewidth',1.5)
%     hold on
    t_fig = figure('Visible','off');
    plot(interp_struct(i).time, interp_struct(i).fluo,'-o', 'Linewidth',1.5)
    title(['Trace ' num2str(i)])
    grid on
%     pause(1.5)
%     hold off
    saveas(t_fig, [outpath,'/traces/', 'trace_' num2str(i) '.png'],'png');
end
%% Run Quality Checks by data set and AP
ap_vec = unique([interp_struct.AP]);
ap_vec = ap_vec(ap_vec > 30);
ap_vec = ap_vec(ap_vec < 50);
% Set size of groups (# AP regions)
grp_size = 4;
ap_grp_ind = floor(ap_vec / grp_size);
ap_grps = unique(ap_grp_ind);
ap_names = {};
for ap_grp = ap_grps
    ap_names = {ap_names{:} [num2str(min(ap_vec(ap_grp_ind==ap_grp))) '-' num2str(max(ap_vec(ap_grp_ind==ap_grp)))]};
end

sets = unique({interp_struct.set});
n_sets = length(sets);
n_grps = length(ap_grps);
maxpause = 0;
for i = 1:length(interp_struct)
    pause_list = [];
    fluo = interp_struct(i).fluo;
    FluoBin = fluo==0;
    ct = 0;
    pct = 0;
    for j = 1:length(fluo)
        ct = ct*FluoBin(j) + FluoBin(j);
        maxpause = max(maxpause, ct);
        if pct > 0 && ct == 0
            pause_list = [pause_list pct];
        end
        pct = ct;
    end
    interp_struct(i).pauses = pause_list;
end
colormap('jet');
cm = colormap;
increment = floor(60 / (n_grps *n_sets));
FluoBins = 1:20:2000;
PauseBins = 1:maxpause;
figure(1);
for i = 1:n_grps
    ap_struct = interp_struct(ismember([interp_struct.AP],ap_vec(ap_grp_ind==ap_grps(i))));
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
        axis([0,1500,0 ,1])
    end
end
hold off
% saveas(figure(1), [outpath, 'fluo_his.png'],'png');
%%
figure(2)
for i = 1:n_grps
    ap_struct = interp_struct(ismember([interp_struct.AP],ap_vec(ap_grp_ind==ap_grps(i))));
    for j = 1:n_sets
        p_list = [];
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                p_list = [p_list ap_struct(a).pauses];
            end
        end
        if isempty(p_list)
            continue
        end
        pause_ct = histc(p_list, PauseBins);
        
        subplot(n_sets,n_grps, (j-1)*n_grps + i)
        bar(PauseBins, pause_ct / max(pause_ct), 'FaceColor',cm(increment*((j-1)*n_grps + i),:),...
            'EdgeColor',cm(increment*((j-1)*n_grps + i),:));
        set(gca,'fontsize',4)
        title(['Pauses, AP: ' ap_names{i}]); %' Set:' sets{j}])
        axis([0,50,0 ,1])
    end
end
saveas(figure(2), [outpath, 'pause_his.png'],'png');
%%
set_titles = {};
for i = 1:length(sets)
    start = strfind(sets{i},'\2017') + 12;
    stop = strfind(sets{i},'\Compiled') - 1;
    string = sets{i};
    set_titles = {set_titles{:}  strrep(string(start:stop),'_',' ')};        
end
ap_range = 39:45;
% set_titles = {'150_1', '200', '150_2', '250 (suspect set)'};
figure(3)
increment = floor(60 / n_sets);
hold on
for j = 1:n_sets
    f_list = [];
    for a = 1:length(interp_struct)
        if strcmp(interp_struct(a).set,sets{j}) && ismember(interp_struct(a).AP,ap_range)
            f_list = [f_list interp_struct(a).fluo];
        end
    end
    if isempty(f_list)
        continue
    end
    ap_ct = histc(f_list, FluoBins);
    plot(FluoBins, cumsum(ap_ct) / sum(ap_ct), 'Color',cm(increment*j,:),'LineWidth', 2);
    
end 
title(['Cumulative PDF (AP: ' num2str(min(ap_range)) '-' num2str(max(ap_range)) ')']); 
axis([0,1500,0,1])
grid on
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(figure(3), [outpath, 'cum_his.png'],'png');

%% Plot Fluo Percentiles By AP and Data Set
ap_vec = unique([interp_struct.AP]);
set_titles = {};
for i = 1:length(sets)
    start = strfind(sets{i},'\2017') + 12;
    stop = strfind(sets{i},'\Compiled') - 1;
    string = sets{i};
    set_titles = {set_titles{:}  strrep(string(start:stop),'_',' ')};        
end

med_fig = figure('Position',[0 0 1024 512]);
increment = floor(60 / n_sets);
medians = nan(n_sets,length(ap_vec));
nTile25 = nan(n_sets,length(ap_vec));
nTile75 = nan(n_sets,length(ap_vec));
% hold on
for j = 1:n_sets
    for k = 1:length(ap_vec)
        f_list = [];
        ap_grp = ap_vec(max(1,k-1):min(length(ap_vec),k+1));
        ap_struct = interp_struct(ismember([interp_struct.AP],ap_grp));
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
ap_vec = unique([interp_struct.AP]);
flank_anterior = 33:38;
flank_posterior = 46:51;
stripe2 = 39:45;
grp_names = {'Anterior Flank', 'Eve2 Stripe', 'Posterior Flank'};
grp_indices = {flank_anterior, stripe2, flank_posterior};
increment = floor(60 /n_sets);
tr_gross = 0;

for k = 1:length(grp_indices)
    fig =  figure('Position',[0 0 1024 512]);
    hold on
    ap_grp = grp_indices{k};
    ap_struct = interp_struct(ismember([interp_struct.AP],ap_grp));
    for j = 1:n_sets
        f_mat = zeros(1,length(time_grid));
        t_list = [];
        tr_ct = zeros(1,length(time_grid));
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                tr_gross = tr_gross + 1;
                f_mat(ismember(time_grid,ap_struct(a).time)) = ...
                    f_mat(ismember(time_grid,ap_struct(a).time)) + ap_struct(a).fluo;                 
                tr_ct = tr_ct + 1*ismember(time_grid,ap_struct(a).time);
            end
        end
        plot(time_grid, f_mat./ tr_ct, 'Color', cm(j*increment,:),'Linewidth',2);
    end
    title(['Mean Fluorescence by Data Set: ' grp_names{k}]) %' (AP: ' num2str(ap_grp) ')'])
    xlabel('seconds');
    ylabel('AU');
    grid on;
    legend(set_titles{:}, 'Location','northeast')
    saveas(fig, [outpath, 'mean_temp_fluo_' grp_names{k} '.png'],'png');
end


%% Histogram of Data Points by AP
figure(6);
subplot(2,1,1);
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

saveas(figure(6), [outpath, 'data_histograms.png'],'png');

% % Truncate Traces with Unreasonably Long Pauses
% % Set reference p_off
% p_off_ref = .9;
% p_cp = 1;
% power = 0;
% while p_cp > .005
%     power = power + 1;
%     p_cp = p_cp * p_off_ref;
% end
% 
% cuts = 0;
% rm_ids = [];
% for i = 1:length(interp_struct)
%     fluo = interp_struct(i).fluo;
%     time = interp_struct(i).time;
%     FluoBin = fluo==0;
%     ct = 0;
%     for j = 1:length(fluo)
%         ct = ct*FluoBin(j) + FluoBin(j);
%         if ct > power
%             fluo_cut = fluo(1:j-ct);
%             time_cut = time(1:j-ct);
%             cuts = cuts + 1;
%             if length(fluo_cut) >= min_len
%                 interp_struct(i).fluo = fluo_cut;
%                 interp_struct(i).time = time_cut;
%                 interp_struct(i).N = length(fluo_cut);
%             else
%                 rm_ids = [rm_ids i];
%             end
%             break
%         end
%     end
% end
% if ~isempty(rm_ids)
%     interp_struct = interp_struct(~ismember(1:length(interp_struct),rm_ids));
% end
% display([num2str(cuts) ' traces truncated due to infeasibly long pauses']);
% 
%save([outpath outName], 'interp_struct');
% 
