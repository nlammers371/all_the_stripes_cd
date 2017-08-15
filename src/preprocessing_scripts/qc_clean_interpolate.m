%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
% folder_path = 'C:\Users\Nicholas\Dropbox (Garcia Lab)\mHMM\orig\';
folder_path = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes';
%Assign Project Identifier
project = 'mHMMeve2_orig_inf';

figpath = ['../../fig/' project '/'];
datapath = ['../../dat/' project '/'];
ap_pos_path = [figpath 'ap_positioning/'];
if exist(ap_pos_path) ~= 7
    mkdir(ap_pos_path)
end  

% Keyword to ensure only sets from current project are pulled
keyword = '20sec';
%vector of data set numbers to include
include_vec = [9,10,19,20,21,22,23,24,26];

%-------------------------Set Summary Parameters--------------------------%
if exist(datapath) ~= 7
    mkdir(datapath);
end
if exist(figpath) ~= 7
    mkdir(figpath);
end
%Pull CompiledParticle Sets from all relevant project folders
dir_struct = struct;
i_pass = 1;
dirinfo = dir(folder_path);
%remove non-directories
dirinfo(~[dirinfo.isdir]) = [];  
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    % skip files lacking project keyword or containing names of skip sets
    if isempty(strfind(thisdir,keyword)) 
        continue
    end
    set_num_start_ind = strfind(thisdir,'_');
    set_num_start_ind = set_num_start_ind(end);
    set_num = str2num(thisdir(set_num_start_ind+1:end));
    
    if sum(set_num==include_vec) ~= 1 
        continue
    end
    
    particle_struct = dir(fullfile(folder_path,thisdir, 'CompiledParticles*'));    
    dir_struct(i_pass).particles = particle_struct;
    dir_struct(i_pass).folder = thisdir;
    i_pass = i_pass + 1;
end
%Save Filepaths of all CompiledParticle Sets
filenames = {};
APnames = {};
for i = 1:length(dir_struct)
    particle_struct = dir_struct(i).particles;
    if length(particle_struct) > 1
        warning(['Multiple Compile Particles Sets Detected for set ' num2str(include_vec(i))])
    end
    filenames = [filenames {[folder_path dir_struct(i).folder '\' particle_struct.name]}];
    APnames = [APnames {[folder_path dir_struct(i).folder '\APDetection.mat']}];
end
%Data structure to store extracted trace sets
trace_struct = struct;
j_pass = 1;
for i = 1:length(filenames)
    %AP Info
    load(APnames{i})
    %Angle between the x-axis and the AP-axis
    APAngle = round(atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)))*360 / (2*pi));    
    %Correction for if APAngle is in quadrants II or III
    if coordPZoom(1)-coordAZoom(1) < 0
        APAngle = APAngle + pi;
    end
    raw_data = load([filenames{i}]);
    time_raw = raw_data.ElapsedTime*60;
    particles = raw_data.CompiledParticles;
    traces_raw = raw_data.AllTracesVector;
    frames_raw = 1:length(time_raw);
    %Get frame that marks start of nc14
    start_frame = raw_data.nc14;
    %Filter trace, frame,and time arrays
    traces = traces_raw(start_frame:end,:);
    time = time_raw(start_frame:end);
    time = time - min(time);
    frames = frames_raw(start_frame:end);
    fn = filenames{i};
    underscores = strfind(fn, '_');
    cp = strfind(fn, '\Compiled');
    %Iterate through traces
    for j = 1:size(traces,2)
        raw_trace = traces(:,j);     
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');
        trace_full = raw_trace(trace_start:trace_stop)';
        trunc_trace = raw_trace(~isnan(raw_trace));
        %Remove dots
        if length(trunc_trace) < 2
            continue
        end        
        %remove NaN values
        trunc_time = time(~isnan(raw_trace));
        time_full = time(trace_start:trace_stop);
        trunc_frames = frames(~isnan(raw_trace));
        %Why the hell is there an APPos field with negative values?
        ap_positions_raw = particles(j).APpos;
        fov_xPos_raw = particles(j).xPos;
        fov_yPos_raw = particles(j).yPos;
        pt_frames = particles(j).Frame;
        ap_positions = ap_positions_raw(ismember(pt_frames,trunc_frames));
        xPos = fov_xPos_raw(ismember(pt_frames,trunc_frames));
        yPos = fov_yPos_raw(ismember(pt_frames,trunc_frames));
        trace_struct(j_pass).APAngle = APAngle;
        trace_struct(j_pass).xPos = xPos;        
        trace_struct(j_pass).yPos = yPos;        
        trace_struct(j_pass).ap_vector = ap_positions;        
        trace_struct(j_pass).fluo = trunc_trace';
        trace_struct(j_pass).time = trunc_time;
        trace_struct(j_pass).fluo_full = trace_full;
        trace_struct(j_pass).time_full = time_full;
        trace_struct(j_pass).setID = str2num(fn(underscores(end)+1:cp-1));
        trace_struct(j_pass).set = filenames{i};
        trace_struct(j_pass).FluoError = particles(j).FluoError;
        j_pass = j_pass + 1;
    end
end

%% Map total fluorecence to AP position for each set Using maximal AP precision (.0001)
%By Default Let's Assume that we wish to take full nc14
start_time = 0*60;
stop_time = 60*60;
%Create Arrays to store total fluorecence and # data points
%Careful with rounding errors here...
ap_index = round(min([trace_struct.ap_vector]),4):.0001:round(max([trace_struct.ap_vector]),4);
ap_index = round(ap_index,4);
ap_fluo_levels = zeros(length(include_vec),length(ap_index));
ap_tp_counts = zeros(length(include_vec),length(ap_index));

%Iterate throughsets
for i = 1:length(include_vec)
    d_set = include_vec(i);
    %Filter for relevant traces
    traces = trace_struct([trace_struct.setID]==d_set);
    %For each time poin in each trace, assign F(t) to AP(t)
    for j = 1:length(traces)
        ap_path = traces(j).ap_vector;
        fluo = traces(j).fluo;
        time = traces(j).time;
        filter_vec = (time >= start_time).*(time < stop_time);
        fluo_trunc = fluo(1==filter_vec); 
        ap_trunc = ap_path(1==filter_vec); 
        for k = 1:length(fluo_trunc)
            ap = round(ap_trunc(k),4);
            %Add fluroescence value to appropriate bin
            ap_filter = ap_index==ap;
            %Make sure that time point is assigned to 1 and only 1 AP
            %position
            if sum(ap_filter) ~= 1
                error('Error in AP Index Assignment');
            end
            ap_fluo_levels(i,ap_filter) = ap_fluo_levels(i,ap_filter) + fluo_trunc(k);
            %Add tp count to appropriate AP bin
            ap_tp_counts(i,ap_filter) = ap_tp_counts(i,ap_filter) + 1;%/length(fluo_trunc);            
        end
    end
end
%% Find Stripe Centers and Save
%Set percent inclusion threshold. Option to test multiple inclusion
%thresholds
pct = .5;
%Set bin size (in frac AP)
ap_bin_size = .02;
%Structure to store results for each data set
stripe_pos_struct = struct;
%Determine stripe center for each dataset
for i = 1:(length(include_vec)+1)
    if i > length(include_vec)
        fluo_vec = sum(ap_fluo_levels);
    else
        fluo_vec = ap_fluo_levels(i,:);
    end
    start = find(fluo_vec>0,1);
    stop = find(fluo_vec,1,'last');
    total_fluo = sum(fluo_vec);
    position_vec = start:stop;
    radius_array = zeros(length(pct),length(position_vec));
    %Find size of stripe for each position
    i_pass = 1;
    for position = start:stop
        fluo_total = fluo_vec(position);
        threshold = pct*total_fluo;
        thresh = threshold;
        %adjust for positions with edges beyond boundary of FOV
        factor = 1;
        radius = 0;
        while fluo_total < thresh
            radius = radius + 1;
            if radius > stop-start
                error('infeasible radius reached')
            end
            factor = (min(stop,position+radius)-max(start,position-radius))/(2*radius+1);
            thresh = factor*threshold;
            fluo_total = sum(fluo_vec(max(start,position-radius):min(stop,position+radius)));
        end
        radius_array(j,i_pass) = radius;
        i_pass = i_pass + 1;
    end
    [m, p] = min(radius_array(j,:));
    best_pos = ap_index(start + p - 1);
    best_rad = m;
    stripe_pos_struct(i).radius_array = radius_array;
    stripe_pos_struct(i).best_pos = best_pos;
    stripe_pos_struct(i).best_rad = best_rad;
    stripe_pos_struct(i).pct = pct;
    if i > length(include_vec)
        stripe_pos_struct(i).setID = 0;
    else
        stripe_pos_struct(i).setID = include_vec(i);
    end
end
% Add Stripe info to trace_struct

for i = 1:length(trace_struct)
    setID = trace_struct(i).setID;
    stripe_info = stripe_pos_struct(include_vec==setID);
    r = stripe_info.best_rad/10000;
    p = stripe_info.best_pos;
    trace_struct(i).stripe_pct = stripe_info.pct;
    trace_struct(i).stripe_rad = r;
    trace_struct(i).stripe_pos = p;
    trace_struct(i).ap_bin_size = ap_bin_size;
    mean_ap = mean(trace_struct(i).ap_vector);
    trace_struct(i).mean_ap = mean_ap;
    ap_diff = mean_ap-trace_struct(i).stripe_pos + ap_bin_size/2;
    binID = int8(sign(ap_diff)*ceil(abs(ap_diff/ap_bin_size)));
    binID(binID>0) = binID(binID>0) - 1; 
    trace_struct(i).APbinID = binID;
end
save([datapath 'raw_traces_' project '.mat'],'trace_struct') 
%% Plot Mean and Cumulative Fluorescence by Set
%Make Titles for Plots (This will lonly work for eve2 format set titles
%currently)
  
sets = unique({trace_struct.set});
n_sets = length(sets);
set_titles = {};
for i = 1:length(sets)
    start = strfind(sets{i},'sec_') + 4;
    stop = strfind(sets{i},'\Compiled') - 1;
    string = sets{i};
    set_titles = {set_titles{:}  string(start:stop)};        
end
%Set dimensions for figs (Need to automate this)
dims = ceil(sqrt(length(include_vec)));
% Set granularity (.01 = 1% AP precision)
precision = 0.01;

%Make Color Palettes for use in figures
colormap('jet')
cm = colormap;
increment = floor(60 / n_sets);
%Array to store color mappings
set_colors = zeros(n_sets, 3);
for i = 1:n_sets
    set_colors(i,:) = cm(1+(i-1)*increment,:);
end

%Generate "coarse" AP vectors to aggregate fine-grained results
%Have to be careful with floating point errors here
coarse_ap_index = round(floor(round(ap_index,4)/precision+10e-6)*precision,2);
coarse_ap = round(unique(coarse_ap_index),2);
coarse_f_avg = zeros(size(ap_fluo_levels,1),length(coarse_ap));
coarse_f_cum = zeros(size(ap_fluo_levels,1),length(coarse_ap));
%Aggregate AP Fluo levels and TP Counts
for i = 1:length(coarse_ap)
   ap = coarse_ap(i);
   coarse_filter = coarse_ap_index==ap;
   if sum(coarse_filter) > precision / .0001 || sum(coarse_filter) < 1
       error('Problem with AP Coarse Graining');
   end
   coarse_f_avg(:,i) = sum(ap_fluo_levels(:,coarse_ap_index==ap),2) ./ sum(ap_tp_counts(:,coarse_ap_index==ap),2);
   coarse_f_cum(:,i) = sum(ap_fluo_levels(:,coarse_ap_index==ap),2);
end

%-------------------------AP Averages-------------------------------------%
mean_fluo_fig = figure('Position',[0 0 1024 1024]);
% Struct to store hist infor for subsequent use
max_mean = max(max(coarse_f_avg));
for j = 1:n_sets        
    subplot(dims,dims,j)
    hold on
    bar(coarse_ap, nanmean(coarse_f_avg)  , 'FaceColor','black',...
        'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);
    bar(coarse_ap, coarse_f_avg(j,:)  , 'FaceColor',set_colors(j,:),...
        'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);
    set(gca,'fontsize',4)
    title(['Mean Fluorescence per Time Step NC 14, Set: ' set_titles{j}]); %' Set:' sets{j}])
    axis([min(coarse_ap),max(coarse_ap),0 ,max_mean])    
    grid on
end
saveas(mean_fluo_fig, [ap_pos_path, 'mean_fluo_ap.png'],'png');

% Integrated Fluorescence With Stripe Centers
cumulative_fluo_fig = figure('Position',[0 0 1024 1024]);
max_cum = 1.1*max(max(coarse_f_cum));
for j = 1:n_sets        
    c = stripe_pos_struct(j).best_pos;
    r = stripe_pos_struct(j).best_rad/10000;
    ac = stripe_pos_struct(end).best_pos;
    subplot(dims,dims, j)
    hold on     
    %Plot average stripe center
    ap = plot([ac,ac],[0,max_cum], 'Linewidth', 1);
    set(ap, 'color', [0.5 0.5 0.5])
    %Plot stripe center for current set
    plot([c,c],[0,max_cum], 'black', 'Linewidth', 1)
    %Make patch indicating stripe size
    p = patch([c-r c+r c+r c-r],[0 0 max_cum max_cum],'black');
    set(p, 'FaceColor',set_colors(j,:),'FaceAlpha', .3)
    bar(coarse_ap, nanmean(coarse_f_cum)  , 'FaceColor','black',...
        'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);        
    bar(coarse_ap, coarse_f_cum(j,:)  , 'FaceColor',set_colors(j,:),...
        'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);        
    set(gca,'fontsize',4)
    title(['Cumulative Fluorescence w/ Stripe (Shaded) in NC 14, Set: ' set_titles{j}],'Fontsize',6); %' Set:' sets{j}])
    axis([min(coarse_ap),max(coarse_ap),0 ,max_cum])    
    grid on
end
saveas(cumulative_fluo_fig, [ap_pos_path, 'cumulative_fluo_ap.png'],'png');
% Temporal Lensing
colormap('winter')
cm = colormap;
ap_fluo_levels = zeros(length(include_vec),length(ap_index));
ap_tp_counts = zeros(length(include_vec),length(ap_index));
ap_nc_counts = zeros(length(include_vec),length(ap_index));

ap_course_index = min(floor(100*ap_index)):max(floor(100*ap_index));
ap_course_index = ap_course_index/100;

start_time = 0*60;
stop_time = 60*50;
window_size = 10*60;
windows = start_time:window_size:stop_time;
n_windows = length(windows)-1;
inc = floor(60/n_windows);

refinement_path = [ap_pos_path '/' 'stripe_refinement/'];
if exist(refinement_path) ~= 7
    mkdir(refinement_path);
end
for i = 1:length(include_vec)
    d_set = include_vec(i);
    traces = trace_struct([trace_struct.setID]==d_set);
    ap_fluo_levels = zeros(n_windows,length(ap_course_index));
    ap_tp_counts = zeros(n_windows,length(ap_course_index));
    ap_nc_counts = zeros(n_windows,length(ap_course_index));
    for j = 1:n_windows
        ss = start_time + (j-1)*window_size;
        e = start_time + j*window_size;
        for k = 1:length(traces)
            ap_path = floor(100*traces(k).ap_vector)/100;
            fluo = traces(k).fluo;
            time = traces(k).time;
            fluo_trunc = fluo(1==(time >= ss).*(time < e)); 
            ap_trunc = ap_path(1==(time >= ss).*(time < e)); 
            for h = 1:length(fluo_trunc)
                ap = ap_trunc(h);
                %Add fluroescence value to appropriate bin
                ap_fluo_levels(j,ap_course_index==ap) = ap_fluo_levels(j,ap_course_index==ap) + fluo_trunc(h);
            end
        end
    end
    ap_fluo_levels = ap_fluo_levels ./ repmat(max(ap_fluo_levels')',1,size(ap_fluo_levels,2));
    ls = {};
    close all
    refine_fig = figure('Position',[0 0 512 1024]);
    for m = 1:size(ap_fluo_levels,1)        
        ind = n_windows - m + 1;
        subplot(n_windows,1,ind);
        bar(ap_course_index,ap_fluo_levels(ind,:),'FaceColor',cm(1+(m-1)*inc,:),'FaceAlpha',.3,'BarWidth', 1)
        title(['Normalized Cumulative Fluorescence Set: ' set_titles{i} ', Minutes ' ...
            num2str(windows(ind)/60) ' to ' num2str(windows(ind+1)/60)]);
        axis([min(ap_course_index) max(ap_course_index) 0 1])
        grid on
    end
    saveas(refine_fig,[refinement_path 'stripe_refinement_set_' set_titles{i} '.png'],'png')
end
%%
%Mean Fluo Over Time
%Structure to store total fluo and fluo counts for each AP and Set
mean_fluo_struct = struct;
max_time = floor(max([trace_struct.time]/60));
min_ap = 34;
max_ap = 45;
bin_size = 2;
all_ap = (min_ap+.5):2:(max_ap-.5);
aggregate_ap_counts = zeros(max_time+1,length(all_ap));
aggregate_ap_fluo = zeros(max_time+1,length(all_ap));
for k = 1:length(include_vec)
    set_traces = trace_struct([trace_struct.setID]==include_vec(k));    
    ap_counts = zeros(max_time+1,length(all_ap));
    ap_fluo = zeros(max_time+1,length(all_ap));
    for j = 1:length(set_traces)       
        t = floor(set_traces(j).time/60);
        f = set_traces(j).fluo;
        ap = floor(100*set_traces(j).ap_vector);
        for i = 1:length(f)
            if ismember(ap(i),floor(all_ap)) || ismember(ap(i),ceil(all_ap)) 
                ap_counts(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap)))...
                    = ap_counts(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap))) + 1;
                ap_fluo(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap))) ...
                    = ap_fluo(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap))) + f(i);
            end
        end
    end
    aggregate_ap_counts = aggregate_ap_counts + ap_counts;
    aggregate_ap_fluo = aggregate_ap_fluo + ap_fluo;
    mean_fluo_struct(k).ap_vec = all_ap;
    mean_fluo_struct(k).ap_counts = ap_counts;
    mean_fluo_struct(k).ap_fluo = ap_fluo;
    mean_fluo_struct(k).set_id = include_vec(k);
end
% Get average overal behavior
aggregate_fluo_mean = sum(aggregate_ap_fluo,2) ./ sum(aggregate_ap_counts,2);

% Make Fig
close all
temporal_fluo_fig = figure('Position',[0 0 2048 2048]);
n_sets = length(include_vec);

increment = floor(60/length(all_ap));
index_vector = 1 + increment*(all_ap - min(all_ap))/2;
colormap('winter');
cm = colormap;

legend_string = {};

for i = 1:n_sets
    subplot(dims,dims,i);
    hold on
    mean_fluo = mean_fluo_struct(i).ap_fluo ./ mean_fluo_struct(i).ap_counts;
    for j = 1:size(mean_fluo,2)
        plot(mean_fluo(:,j),'Color',[cm(index_vector(all_ap(j)==all_ap),:) .5],'Linewidth',1.5);            
    end
    plot(aggregate_fluo_mean, '-', 'Color', 'black','Linewidth',1.5)
    axis([0 50 0 1000])
    grid on
    xlabel('Minutes into nc14');
    ylabel('Mean Fluorescence (AU)');
    title(['Mean Fluorescence Over Time by AP: Set ' num2str(include_vec(i))],'FontSize',8);
    set(gca,'FontSize',6)
end

for ap = all_ap
    ap_string = ['AP ' num2str(ap)];    
    legend_string = {legend_string{:} ap_string};
end
legend(legend_string{:});

saveas(temporal_fluo_fig, [ap_pos_path, 'mean_fluo_temp.png'],'png');
%% Multi Fluo Hist Plots
close all

%Set size of fluo bins
granularity = 20;
max_fluo = ceil(max([trace_struct.fluo])/granularity)*granularity;
FluoBins = 1:granularity:max_fluo;
fluopath = [figpath '/fluo_statistics/'];
if exist(fluopath) ~= 7
    mkdir(fluopath)
end
fluo_his_fig = figure('Position',[0 0 1024 1024]);
% Struct to store hist infor for subsequent use
hist_info = struct;

for j = 1:n_sets
    bin_struct = trace_struct([trace_struct.binID]==0);
    f_list = [];
    for a = 1:length(bin_struct)
        if bin_struct(a).setID==include_vec(j)
            fluo = bin_struct(a).fluo;
            f_list = [f_list fluo(fluo>0)];
        end
    end
    if isempty(f_list)
        continue
    end
    ap_ct = histc(f_list, FluoBins);        
    subplot(dims,dims, j)
    b = bar(FluoBins, ap_ct / max(ap_ct), 'FaceColor',set_colors(j,:),...
        'EdgeColor',set_colors(j,:),'BarWidth', 1);
    set(gca,'fontsize',4)
    title(['Fluo Distribution in Stripe Center, Set: ' set_titles{j}]); 
    axis([0,max_fluo,0 ,1])    
    grid on
end

saveas(fluo_his_fig, [fluopath 'set_stripe_fluo.png'],'png');
hold off

%% Cumulative Fluorescence within Eve Stripe 2 Region
close all
%Make Strings for legen entries
cf_fig = figure;
ptile_list = zeros(1,n_sets);
% Set scale for x axis. Max of selected percentile across sets will be used
ptile = 97;
hold on
for j = 1:n_sets
    f_list = [];
    for a = 1:length(trace_struct)
        if trace_struct(a).setID == include_vec(j) && trace_struct(a).binID==0
            f_list = [f_list trace_struct(a).fluo];
        end
    end
    ptile_list(j) = prctile(f_list,ptile); 
    if isempty(f_list)
        continue
    end
    ap_ct = histc(f_list, FluoBins);
    plot(FluoBins, cumsum(ap_ct) / sum(ap_ct), 'Color',set_colors(j,:),'LineWidth', 2);
end 
title('Cumulative PDF for Stripe Centers');
axis([0,max(ptile_list),0,1])
grid on
xlabel('AU')
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(cf_fig, [fluopath, '/set_cum_fluo.png'],'png');
