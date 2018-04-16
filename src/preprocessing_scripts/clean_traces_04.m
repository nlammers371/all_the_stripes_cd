% Script to clean and interpolate fluorescence trace
% Also clean up nucleus data-set
clear 
close all
%%%Set Cleaning and Summary Parameters
min_dp = 15; % minimum # dp acceptable for an inference trace
Tres_interp = 20;
InterpGrid = 0:Tres_interp:60*50;
FOV_edge_padding = 20; % pixels
xDim = 1024;
%------------------------Import Raw Trace Set------------------------%
%ID's of sets to includeg
project = 'eve7stripes_inf_2018_03_27_final';
print_traces = 0; %Output PNG files for traces?

%---------------------------Set Paths-------------------------------------%
TracePath = ['../../dat/' project '/' 'raw_traces_' project '.mat'];
NucleusPath = ['../../dat/' project '/' 'ellipse_info_' project '.mat'];
OutPath = ['../../dat/' project '/'];
FigPath = ['../../fig/experimental_system/' project '/preprocessing'];
TraceSavePath = [FigPath '/traces/'];
mkdir(TraceSavePath);
mkdir(OutPath);

% Save Names
CleanTraceName = ['inference_traces_' project '_dT' num2str(Tres_interp) '.mat'];
CleanNucleusName = ['inference_nuclei_' project '_dT' num2str(Tres_interp) '.mat'];

%----------------Load Traces and Perform First Pass QC--------------------%
%Load raw traces (saved in struct titled "trace_struct")
load(TracePath);
load(NucleusPath);
set_index = unique([trace_struct.setID]);
%%% Cleaning Parameters
big_jump1 = prctile([trace_struct.fluo],99);
jump_threshold1 = big_jump1/1.5; % this should be a conservative threshold for single time step increase
index_vec = 1:length(trace_struct); % convenience ref vector
field_names = fieldnames(trace_struct);
%%% remove small set of duplicates
rm_indices = [];
particle_vec = [trace_struct.ParticleID];
for p = 1:length(particle_vec)
    pID = trace_struct(p).ParticleID;
    matches = find(particle_vec == pID);
    if length(matches) > 1
        rm_indices = [rm_indices matches(2:end)];
    end
end
trace_struct = trace_struct(~ismember(1:length(trace_struct),rm_indices));

trace_struct_final = []; % Structure to store traces after first round of cleaning
jump_ct = 0;
blip_ct1 = 0;
for i = 1:length(trace_struct) 
    temp = trace_struct(i);
    if temp.x_flip_flag        
        temp.xPos = xDim - temp.xPos + 1;
    end    
    trace1 = temp.fluo; %Load full trace, including intervening NaN's
    time = temp.time;      
    quality_flag = 1;
    if sum(~isnan(trace1)) < min_dp
        quality_flag = 0;
    end
    %Null assumption is that all clusters of 6 or more NaNs are 0s. Single
    %,double, or triple NaNs are assumed to have been missed nonzero dps
    trace1_nans = isnan(trace1);      
    %Look for clusters of 6 or more NaNs
    kernel = [1,1,1,1,1];
    tn_conv = conv(kernel,trace1_nans);
    tn_conv = tn_conv(3:end-2);
    z_ids = find(tn_conv==5);
    z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
    trace1(z_ids) = 0; % set clusters to zeros
    trace1(trace1<0) = 0; % deal with negative values
    % find single dp "blips"
    tr_dd1 = abs([0 diff(diff(trace1)) 0]);  % 1 slice
    trace1(tr_dd1>2*jump_threshold1) = NaN;
    blip_ct1 = blip_ct1 + sum(tr_dd1>2*jump_threshold1);
    % interpolate remaining NaNs    
    query_points1 = time(isnan(trace1));%InterpGrid((InterpGrid>=min(time))&(InterpGrid<=max(time)));
    interp_t1 = time(~isnan(trace1));
    interp_f1 = trace1(~isnan(trace1));
    new_f1 = interp1(interp_t1,interp_f1,query_points1);  
    trace1(ismember(time,query_points1)) = new_f1;    
    %%% flag traces with unreasonably large rises or falls    
    tr_d1 = diff(trace1);
    if max(abs(tr_d1)) >= jump_threshold1 % || max(abs(tr_d3)) >= jump_threshold3
        jump_ct = jump_ct + 1;
        quality_flag = 0;
    end    
    % Interpolate to standardize spacing
    t_start = InterpGrid(find(InterpGrid>=time(1),1));
    t_stop = InterpGrid(find(InterpGrid<=time(end),1,'last'));
    time_interp = t_start:Tres_interp:t_stop;
    trace1_interp = interp1(time,trace1,time_interp);
    interp_fields = {'xPos','yPos','ap_vector','stripe_id_vec'};
    % interpolate other vector fields
    for j = 1:length(interp_fields)
        int_vec = temp.(interp_fields{j});
        int_time = temp.time;
        cp_frames = temp.cp_frames;
        all_frames = temp.all_frames;
        int_time = int_time(ismember(all_frames,cp_frames));
        temp.([interp_fields{j} '_interp']) = interp1(int_time,int_vec,time_interp);
    end     
    temp.stripe_id_vec_interp = round(3*temp.stripe_id_vec_interp)/3; 
    temp.fluo_interp = trace1_interp;
    temp.time_interp = time_interp;
    temp.inference_flag = quality_flag;
    trace_struct_final = [trace_struct_final temp];    
end

% Find last obs times for each set
last_times = [];
for s = 1:length(set_index)
    setID = set_index(s);
    last_times = [last_times max([trace_struct_final([trace_struct_final.setID]...
        ==setID).time_interp])];
end
time_ceiling = InterpGrid(find(InterpGrid<=min(last_times),1,'last')); 
%%% Add a few useful fields
for i = 1:length(trace_struct_final)
    fluo_interp = trace_struct_final(i).fluo_interp;
    time_interp = trace_struct_final(i).time_interp;    
    setID = trace_struct_final(i).setID;        
    time_full = InterpGrid;
    fluo_full1 = zeros(1,length(time_full));
    fluo_full1(ismember(time_full,time_interp)) = fluo_interp;
    trace_struct_final(i).fluo_full = fluo_full1; %"unclipped" version
    trace_struct_final(i).time_full = time_full;    
    trace_struct_final(i).N = length(fluo_interp);
    trace_struct_final(i).dT = Tres_interp;                
    trace_struct_final(i).alpha_frac = 1302/6444;
    trace_struct_final(i).InterpGrid = InterpGrid;
    trace_struct_final(i).FOV_edge_padding = FOV_edge_padding;    
end

%------------------------- Clean Ellipse Set -----------------------------%

padding = 10; % tolerance for missing frames at end or beginning of nc14 
trace_particle_vec = [trace_struct_final.ParticleID];
nuclei_clean = [];
rm_counts = 0;
set_vec = unique([trace_struct_final.setID]);
for i = 1:length(schnitz_struct)
    temp = schnitz_struct(i);
    if temp.x_flip_flag
        temp.xPos = xDim - temp.xPos + 1;
    end
    setID = temp.setID;
    lt = last_times(set_index==setID);
    nc_times = temp.time;         
    temp.InterpGrid = InterpGrid;
    % interpolate relevant fields
    time = temp.time;
    frames = temp.frames;    
    t_start = InterpGrid(find(InterpGrid>=time(1),1));    
    t_stop = InterpGrid(find(InterpGrid<=time(end),1,'last'));
    time_interp = t_start:Tres_interp:t_stop;
    frames_interp = linspace(min(frames), max(frames), length(time_interp));
    tracking_flags = ismember(floor(frames_interp),frames);
    interp_fields = {'xPos','yPos','time','ap_vector','stripe_id_vec'};
    for j = 1:length(interp_fields)
        int_vec = temp.(interp_fields{j});
        int_time = temp.time;
        interp = interp1(int_time,int_vec,time_interp);
        interp(~tracking_flags) = NaN;
        temp.([interp_fields{j} '_interp']) = interp;
    end    
    temp.stripe_id_vec_interp = round(temp.stripe_id_vec_interp*3)/3;
    FOV_flags = (temp.xPos_interp > FOV_edge_padding) & (temp.xPos_interp < 1024- FOV_edge_padding)...
        &(temp.yPos_interp > FOV_edge_padding) & (temp.yPos_interp < 256 - FOV_edge_padding);    
    FOV_flags_all = zeros(size(InterpGrid));
    FOV_flags_all(ismember(InterpGrid,time_interp)) = 1;
    temp.FOV_flags = FOV_flags_all;
    temp.nc14_flag = time(end) >= time_ceiling && time(1) <= 360;
    temp.FOV_edge_padding = FOV_edge_padding;            
    ParticleID = temp.ParticleID;
    trace_ind = [];
    if ~isnan(ParticleID)
        trace_ind = find(trace_particle_vec==ParticleID);
        if isempty(trace_ind) % catch cases in which particle has been removed for QC reasons
            error('Inconsistent indexing between Trace and Nucelus Structures')
        else
            ncIDCheck = trace_struct_final(trace_ind).ncID;
            if ncIDCheck ~= temp.ncID
                error('Inconsistent Particle and NC Identifiers')
            end
        end
    else
        trace_ind = NaN;
    end
    if length(trace_ind) > 1
        error('Degenerate Identifiers')
    end
    temp.TraceIndex = trace_ind;
    nuclei_clean = [nuclei_clean temp];
end
%%% Now add flag to trace struct indicating which traces correspond to
%%% "good" nuclei
% cross-reference id variables for consistency check
trace_nucleus_vec = [trace_struct_final.ncID];
% trace index var
NucleusTraceIndices = [nuclei_clean.TraceIndex];

for i = 1:length(trace_struct_final)
    pID = trace_struct_final(i).ParticleID;
    ncInd = find(NucleusTraceIndices==i);
    pIDCross = nuclei_clean(ncInd).ParticleID;
    if pIDCross ~= pID
        error('Inconsistent Identifiers')
    end
    trace_struct_final(i).FOV_flags = nuclei_clean(ncInd).FOV_flags;
    trace_struct_final(i).nc14_flag = nuclei_clean(ncInd).nc14_flag;
end

% save
save([OutPath CleanTraceName],'trace_struct_final');
save([OutPath CleanNucleusName],'nuclei_clean');

%%% If desired, save individual trace plots
if print_traces
    for i = 1:30:length(trace_struct_final)
        t_fig = figure('Visible','off');
        cm = jet(64);
        hold on
        plot(trace_struct_final(i).time / 60, trace_struct_final(i).fluo...
                , '-o','Linewidth',1.5,'Color',cm(15,:)/1.5)
        plot(trace_struct_final(i).time_interp / 60, trace_struct_final(i).fluo_interp,...
            '-s', 'Linewidth',1.5,'Color',cm(15,:))
        inf_flag = trace_struct_final(i).inference_flag;
        title(['Original vs. Final: Trace ' num2str(i) ' (inf flag=' num2str(inf_flag) ')'])
        legend('Raw', 'Interpolated');
        xlabel('Minutes');
        grid on        
        saveas(t_fig, [TraceSavePath, 'trace_' num2str(i) '.png'],'png');
    end
end