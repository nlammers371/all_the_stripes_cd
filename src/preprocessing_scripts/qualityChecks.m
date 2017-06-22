%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = '../../../data/';
outpath = '../../processed_data/data_quality/';
outName = 'eveSet_2017_06_15.mat';
if exist(outpath) ~= 7
    mkdir(outpath);
end
files = dir(folder_path);

filenames = {};
for i = 1:length(files)
    if strfind(files(i).name,'.mat') > 0
        filenames = [filenames {files(i).name}];
    end
end

traces_no_nan = struct;
traces_interp = struct;
diff_traces = 0;
%Data structure to store extracted trace sets

trace_struct = struct;
i_iter = 1;
for k = 1:length(filenames)
    raw_data  = load([folder_path filenames{k}]);
    time = raw_data.ElapsedTime*60;
    time = time - time(raw_data.nc14);
    traces = raw_data.AllTracesVector;
    nc14_filter = raw_data.ncFilter(:,end);
    traces_14 = traces(:,nc14_filter == 1);
    %Set Nans to -1
    traces_14(traces_14==0) = -1;
    traces_14(isnan(traces_14)) = 0;
    for i = 1:size(traces_14,2)
        raw_trace = traces_14(:,i);
        start = find(raw_trace,1);
        stop = find(raw_trace,1,'last');
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

%Get 90th percentile for point-to-point deltas
ref_len = prctile(abs(diff([trace_struct.fluo])),90);
adjustments = 0;
for i = 1:length(trace_struct)
    trace = trace_struct(i).fluo;
    tr_d = [0 diff(trace)];
    tr_dd = [0 diff(diff(trace)) 0];
    rm_list = [];
    for j = 1:length(trace)
        if abs(tr_d(j)) > .5 * ref_len && trace(j) == 0
            rm_list = [rm_list j];
            adjustments = adjustments +1;
        elseif abs(tr_dd) > 1.5*ref_len
            rm_list = [rm_list j];
            adjustments = adjustments +1;
        end
    end
    trace_struct(i).fluo = trace(~ismember(1:length(trace), rm_list));
    trace_struct(i).time = trace_struct(i).time(~ismember(1:length(trace),rm_list));
end
        
% Interpolate to Achieve Tractable Memory
%Define Desired res and associated param
mem = 8;
T_elong = 150;
Tres = T_elong / mem;
%Set minimum trace length (in time steps)
min_len = mem;
i_pass = 1;
interp_struct = struct;
for i = 1:length(trace_struct)
    fluo = trace_struct(i).fluo;
    time = trace_struct(i).time;
    
    time_interp= min(time) + Tres*(0:floor((max(time)-min(time))/Tres));
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
%% Run Quality Checks by data set and AP
ap_vec = unique([interp_struct.AP]);
% Set size of groups (# AP regions)
grp_size = 5;
ap_grp_ind = floor(ap_vec / grp_size);
ap_grps = unique(ap_grp_ind);
ap_names = {};
for ap = ap_grps
    ap_names = {ap_names{:} [num2str(min(ap_vec(ap_grp_ind==ap))) '-' num2str(max(ap_vec(ap_grp_ind==ap)))]};
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
FluoBins = 1:5000:450000;
PauseBins = 1:maxpause;
figure(1)
for i = 1:n_grps
    ap_struct = interp_struct(ismember([interp_struct.AP],ap_vec(ap_grp_ind==ap_grps(i))));
    for j = 1:n_sets
        f_list = [];
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                f_list = [f_list ap_struct(a).fluo];
            end
        end
        if isempty(f_list)
            continue
        end
        fluo_ct = histc(f_list, FluoBins);
        
        subplot(n_sets,n_grps, (j-1)*n_grps + i)
        bar(FluoBins, fluo_ct / max(fluo_ct), 'FaceColor',cm(increment*((j-1)*n_grps + i),:),...
            'EdgeColor',cm(increment*((j-1)*n_grps + i),:));
        set(gca,'fontsize',4)
        title(['Fluo, AP: ' ap_names{i}]); %' Set:' sets{j}])
        axis([0,400000,0 ,1])
    end
end
saveas(figure(1), [outpath, 'fluo_his.png'],'png');
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
set_titles = {'150_1', '200', '150_2', '250 (suspect set)'};
figure(3)
increment = floor(60 / n_sets);
hold on
for j = 1:n_sets
    f_list = [];
    for a = 1:length(interp_struct)
        if strcmp(interp_struct(a).set,sets{j})
            f_list = [f_list interp_struct(a).fluo];
        end
    end
    if isempty(f_list)
        continue
    end
    fluo_ct = histc(f_list, FluoBins);
    plot(FluoBins, cumsum(fluo_ct) / sum(fluo_ct), 'Color',cm(increment*j,:),'LineWidth', 2);
    title('Cumulative PDF (All AP Positions)'); 
    axis([0,400000,0,1])
end 
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(figure(3), [outpath, 'cum_his.png'],'png');

%%

% Truncate Traces with Unreasonably Long Pauses
%Set reference p_off
p_off_ref = .9;
p_cp = 1;
power = 0;
while p_cp > .005;
    power = power + 1;
    p_cp = p_cp * p_off_ref;
end

cuts = 0;
rm_ids = [];
for i = 1:length(interp_struct)
    fluo = interp_struct(i).fluo;
    time = interp_struct(i).time;
    FluoBin = fluo==0;
    ct = 0;
    for j = 1:length(fluo)
        ct = ct*FluoBin(j) + FluoBin(j);
        if ct > power
            fluo_cut = fluo(1:j-ct);
            time_cut = time(1:j-ct);
            cuts = cuts + 1;
            if length(fluo_cut) >= min_len
                interp_struct(i).fluo = fluo_cut;
                interp_struct(i).time = time_cut;
                interp_struct(i).N = length(fluo_cut);
            else
                rm_ids = [rm_ids i];
            end
            break
        end
    end
end
if ~isempty(rm_ids)
    interp_struct = interp_struct(~ismember(1:length(interp_struct),rm_ids));
end
display([num2str(cuts) ' traces truncated due to infeasibly long pauses']);

%save([outpath outName], 'interp_struct');