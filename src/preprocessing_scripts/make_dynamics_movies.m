% Script to generate movies illustrating promoter and silencing dynamics
%%% Load Data & Create Write Paths
clear all
close all
%------------------------Set Path Specs, ID Vars------------------------%
FolderPath = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
FISHPath = 'D:\Data\Augusto\LivemRNA\Data\PreProcessedData';
project = 'eve7stripes_inf_2018_02_20'; %Project Identifier
% folders
fig_path = ['../../fig/experimental_system/' project '/movies/'];
data_path = ['../../dat/' project '/']; % data mat directory
% fig subfolders
ap_pos_path = [fig_path 'ap_positioning/'];
fluo_path = [fig_path 'fluo_stats/'];
trace_name = [data_path 'raw_traces_' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_' project]; % names for compiled elipse struct

% make filepaths
mkdir(data_path);
mkdir([ap_pos_path '/stripe_fits']);
mkdir(fluo_path);
% cleaning params
keyword = '_30uW_550V'; % Keyword to ensure only sets from current project are pulled
% show_ap_fit_figs = 0;
snippet_size = 15; % particles within snippet/2+1 are of concern
pre_post_padding = 10; % max mun frames for which nucleus can be MIA at start or end
% store set names
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
cp_filenames = {}; % particles
ellipse_filenames = {}; % ap info
nc_filenames = {}; % nuclei
set_nums = [];
prefix_list = {};
frame_list = {};
for d = 1 : length(dirinfo)
    thisdir = dirinfo(d).name;
    % Skip files lacking project keyword 
    if isempty(strfind(thisdir,keyword)) 
        continue
    end    
    prefix_list = [prefix_list{:} {thisdir}];
%     set_num_start_ind = strfind(thisdir,'_');
%     set_num_start_ind = set_num_start_ind(end);
%     set_num = str2num(thisdir(set_num_start_ind+1:end));        
%     set_nums = [set_nums set_num];
    % append file paths
    cp_name = dir([FolderPath '/' thisdir '/CompiledParticles*']);
    cp_name = cp_name(1).name;
    cp_filenames = [cp_filenames {[thisdir '/' cp_name]}];    
    frame_list = [frame_list{:} {[thisdir '/FrameInfo.mat']}];
    ellipse_filenames = [ellipse_filenames {[thisdir '/Ellipses.mat']}];    
    nc_filenames = [nc_filenames {[thisdir '/' thisdir '_lin.mat']}];           
end
cm = jet(128);
increment = floor(size(cm,1) / 7);
%Array to store color mappings
stripe_colors = zeros(7, 3);
for i = 1:7
    stripe_colors(i,:) = cm(1+(i-1)*increment,:);
end


%Embryo ID and Folder location
for set_num = 1:length(cp_filenames)
    Prefix = prefix_list{set_num};

    %Load the data
    load([FolderPath '/' frame_list{set_num}])
    load([FolderPath '/' cp_filenames{set_num}])
    load([FolderPath '/' nc_filenames{set_num}])
    load([FolderPath '/' ellipse_filenames{set_num}])

    %%% ------ movie parameters
    MaxRadius = 22;
    MaxFluo = 250000;
    n_sets = length(cp_filenames);
    xDim = FrameInfo.PixelsPerLine;
    yDim = FrameInfo.LinesPerFrame;
    [px, py] = meshgrid(1:xDim,1:yDim);
    CurrentNC = 14;

    %%% load inference and inference traces results (needed for spot state info)
    load('D:\Data\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_02_20\raw_traces_eve7stripes_inf_2018_02_20.mat')

    %%% ------------------- Generate Movies ------------------------------- %%%
    setID = set_num;
    last_frame = max([trace_struct([trace_struct.setID]==setID).all_frames]);
    % first_frame = min([trace_struct_filtered([trace_struct_filtered.setID]==setID).all_frames]);
    FrameRange=nc14:last_frame;

    %%% Make Movie with Nucleus Masks Indicating Instantaneous Spot Fluo
    % Set write directory
    fluo_path = [fig_path '/' Prefix '/fluo_masks/'];
    mkdir(fluo_path);
    set_struct = trace_struct([trace_struct.setID]==set_num);
    % make ref vectors
    fluo_vec = [];
    time_vec = [];
    frame_vec = [];
    nucleus_vec = [];
    stripe_id_vec = [];
    xPos_vec = [];
    for j = 1:length(set_struct)
        fluo = set_struct(j).fluo;
        time = set_struct(j).time;
        time = time(~isnan(fluo));
        fluo = fluo(~isnan(fluo));   
        frames = set_struct(j).cp_frames;
        xPos = set_struct(j).xPos;
        time_vec = [time_vec time];
        fluo_vec = [fluo_vec fluo];
        frame_vec = [frame_vec frames];
        nucleus_vec = [nucleus_vec repelem(set_struct(j).Nucleus,length(fluo))];
        stripe_id_vec = [stripe_id_vec repelem(set_struct(j).stripe_id_coarse,length(fluo))];
        xPos_vec = [xPos_vec xPos];
    end
    max_stripe = nanmax(stripe_id_vec);
    min_stripe = nanmin(stripe_id_vec);
    mean_x_min = mean(xPos_vec(stripe_id_vec==min_stripe));
    mean_x_max = mean(xPos_vec(stripe_id_vec==max_stripe));
    if mean_x_max > mean_x_min
        flip_flag = 0;
    else
        flip_flag = 1;
    end
    %Iterate Through Frames
    for CurrentFrame = FrameRange    
        %Highlight Active Regions with Viterbi State    
        %Track pixel assignments
        NucleusFluoMat = NaN(yDim,xDim,3);
        NucleusIDMat = zeros(yDim,xDim);
        NucleusDistMat = ones(yDim,xDim)*MaxRadius;    
        % Loop through ALL nuclei (inactive and active) and assign pixels
        % to a nucleus
        for s = 1:length(schnitzcells)
            try
                CurrentEllipse= schnitzcells(s).cellno(...
                                schnitzcells(s).frames==...
                                CurrentFrame);
                x =  Ellipses{CurrentFrame}(CurrentEllipse,1)+1;
                y =  Ellipses{CurrentFrame}(CurrentEllipse,2)+1;        
                if isempty(x)
                    continue
                end 
                distances = ((px-x).^2 + (py-y).^2).^.5;
                candidate_indices = NucleusDistMat > distances; 
                %Record Fluorescence        
                NucleusIDMat(candidate_indices) = s;
                NucleusDistMat(candidate_indices) = distances(candidate_indices);
            catch
                warning(['Nucleus Indexing Error in Frame ' num2str(CurrentFrame)])
            end
        end
        % reference frame index cell
        extant_nuclei = nucleus_vec(frame_vec==CurrentFrame);
        extant_fluo = fluo_vec(frame_vec==CurrentFrame);
        extant_stripes = stripe_id_vec(frame_vec==CurrentFrame);
        ParticlesToShow = [];
        for i = 1:length(extant_nuclei)
            Nucleus = extant_nuclei(i);     
            fluo = extant_fluo(i);
            stripe_id = extant_stripes(i);
            if isnan(stripe_id)
                stripe_color = [1 1 1]*min(1,max(.4,fluo/MaxFluo));
            else
                stripe_color = stripe_colors(stripe_id,:)*min(1,max(.2,fluo/MaxFluo));
            end
            CurrentEllipse=...
                schnitzcells(Nucleus).cellno(...
                schnitzcells(Nucleus).frames==...
                CurrentFrame);                    
            %Record Fluorescence
            for k = 1:3
                filter = NucleusIDMat==Nucleus;        
                nc_slice = NucleusFluoMat(:,:,k);
                nc_slice(filter) = stripe_color(k);
                NucleusFluoMat(:,:,k) = nc_slice;        
            end
            ParticlesToShow = [ParticlesToShow Nucleus];
    %             display(['Error in frame ',num2str(CurrentFrame)])                                         
        end

        %Now Draw Nucleus Borders for active nuclei
        NucleusBorderMat = zeros(size(NucleusIDMat));
        window = 1;
        for i = ParticlesToShow
            %get coordinates of nucleus patch
            if sum(sum(NucleusIDMat==i)) > 0
                x_vec = reshape(px(NucleusIDMat==i),[],1);
                y_vec = reshape(py(NucleusIDMat==i),[],1);
                for j = 1:length(x_vec)
                    metric = sum(sum(NucleusIDMat(max(1,y_vec(j)-window):min(y_vec(j)+window,yDim),...
                             max(1,x_vec(j)-window):min(x_vec(j) + window,xDim))));
                    if metric~= i*(2*window+1)^2
                        NucleusBorderMat(y_vec(j),x_vec(j)) = 1;
                    end
                end
            end
        end
        %Prevent overlap between fluroescence mask and borders
        max_p = max(NucleusFluoMat,[],3);
        for k = 1:3
            nc_slice = NucleusFluoMat(:,:,k);
%             nc_slice(isnan(max_p)) = .2;
            nc_slice(NucleusBorderMat>0) = 0;
            NucleusFluoMat(:,:,k) = nc_slice;
        end
    %     %Make a maximum projection of the mRNA channel
    %     D=dir([FISHPath,'\',Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z*.tif']);
    %     %Do not load the first and last frame as they are black
    %     ImageTemp=[];
    %     for m=2:(length(D)-1)
    %         ImageTemp(:,:,m-1)=imread([FISHPath,'\',Prefix,filesep,D(m).name]);
    %     end
    %     mRNAImage=max(ImageTemp,[],3);        
    %     %Load the corresponding histone image
    %     HistoneImage=imread([FISHPath,'\',Prefix,filesep,...
    %         Prefix,'-His_',iIndex(CurrentFrame,3),'.tif']);        
    % %     HistoneImage = HistoneImage/20;
    %     %find most likely state
    %     
    %     %Overlay all channels       

        OverlayFig = figure('Visible','off');
        clf
        if flip_flag
            imshow(fliplr(NucleusFluoMat))   
        else
            imshow(NucleusFluoMat)   
        end

        text(20,225,[iIndex(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),2),...
            ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])
        xlim([0,xDim])    
        ylim([0,yDim])
        drawnow

        SavedFlag=0;

%         while ~SavedFlag
%             try
%                 saveas(gcf,[fluo_path '\nc',num2str(CurrentNC),...
%                     '-',num2str(iIndex(CurrentFrame,3)), '.tif']);   
%                 SavedFlag=1;
%             catch
%                 disp('Error saving. Retrying...')
%             end
%         end    
    end
    close all
end