% Script to extract anterior and posterior coordinates for each data set
close all
clear 
%------------------------Set Path Specs, ID Vars------------------------%
FolderPath = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
project = 'eve7stripes_inf_2018_04_28'; %Project Identifier
% folders
fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
out_path = ['../../dat/' project '/'];
% cleaning params
keyword = '_30uW_550V'; % Keyword to ensure only sets from current project are pulled

% store set names
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
ap_filenames = {}; % ap info
set_nums = [];
for d = 1 : length(dirinfo)
    thisdir = dirinfo(d).name;
    % Skip files lacking project keyword 
    if isempty(strfind(thisdir,keyword)) 
        continue
    end    
    % append file paths    
    ap_filenames = [ap_filenames {[thisdir '/APDetection.mat']}];        
end

coord_set = NaN(length(ap_filenames),5);
header = {'set_id', 'a_x','a_y', 'p_x', 'p_y'};
for i = 1:length(ap_filenames) % Loop through filenames    
    % read in raw files
    load([FolderPath ap_filenames{i}]) % AP Info   
    coord_set(i,1) = i;
    coord_set(i,2:3) = CoordA;
    coord_set(i,4:5) = CoordP;
end