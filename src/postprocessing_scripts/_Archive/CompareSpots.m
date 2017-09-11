%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\';
sub_paths = {'Eve2_ML', 'Eve2_orig'};
key_names = {'ML', 'ORIG'};
project = 'mHMMeve2_ml_comparisons_07_12';
outpath = [folder_path '/projects/' project '/' ];
% Keyword to ensure only sets from current project are pulled
keyword = 'eve2_20sec_';
exclude = 'eve2_20sec_5';
if exist(outpath) ~= 7
    mkdir(outpath);
end
%%
meta_struct = struct;
for h = 1:length(key_names)
    dir_struct = struct;
    i_pass = 1;
    dirinfo = dir([folder_path sub_paths{h}]);
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    subdirinfo = cell(length(dirinfo));
    for K = 1 : length(dirinfo)
      thisdir = dirinfo(K).name;
      % skip files lacking prject keyword
      if isempty(strfind(thisdir,keyword)) || ~isempty(strfind(thisdir,exclude))
          continue
      end
      subdir_struct = dir(fullfile(folder_path,sub_paths{h},thisdir, '*.mat'));
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
    spot_names = {};
    for i = 1:length(dir_struct)
        subdir_struct = dir_struct(i).files;
        file_count = 0;
        filenames_old = filenames;
        spot_names_old = spot_names;
        for j = 1:length(subdir_struct)
            if strfind(subdir_struct(j).name,'CompiledParticles') > 0
                filenames = [filenames {[subdir_struct(j).folder '\' subdir_struct(j).name]}];
                file_count = file_count + 1;
            end
            if strfind(subdir_struct(j).name,'Spots') > 0
                spot_names = [spot_names {[subdir_struct(j).folder '\' subdir_struct(j).name]}];
                file_count = file_count + 1;
            end
        end
        if file_count ~= 2
            filenames = filenames_old;
            spot_names = spot_names_old;
        end
    end
    
    %Data structure to store extracted trace sets
    trace_struct = struct;
    spot_struct = struct;
    particle_struct = struct;
    i_iter = 1;
    for k = 1:length(filenames)
        raw_data  = load([filenames{k}]);
        spot_struct = load([spot_names{k}]);
        spot_struct = spot_struct.Spots;
        time = raw_data.ElapsedTime*60;
        traces = raw_data.AllTracesVector;
        %Set Nans to 0, true zeros to -1
        traces(traces==0) = -1000;
        traces(isnan(traces)) = 0;
        particles = raw_data.CompiledParticles;
        fName = filenames{k};
        setID = str2num(fName(strfind(fName,'20sec_')+6:strfind(fName,'\Comp')-1));
        for i = 1:size(traces,2)
            raw_trace = traces(:,i);
            nucleus = particles(i).Nucleus;
            particle = particles(i).OriginalParticle;
            start = find(raw_trace,1);
            stop = find(raw_trace,1,'last');
            raw_trace(raw_trace==0) = nan;
            raw_trace(raw_trace==-1000) = 0;
            trunc_trace = [raw_trace(start:stop)'];
            trunc_time = time(start:stop); 
            [~, apPos] = max(raw_data.APFilter(i,:));
            trace_struct(i_iter).Index = particles(i).Index;
            trace_struct(i_iter).fluo = trunc_trace;
            trace_struct(i_iter).Frame = particles(i).Frame;
            trace_struct(i_iter).time = trunc_time;
            trace_struct(i_iter).setID = setID;
            trace_struct(i_iter).ncID = nucleus;
            trace_struct(i_iter).pID = particle;
            i_iter = i_iter + 1;
        end
        meta_struct(k).(['raw_' key_names{h}]) = trace_struct;
        meta_struct(k).(['spots_' key_names{h}]) = spot_struct;
        meta_struct(k).(['particles_' key_names{h}]) = particles;
    end
end

%%
t_pass = 1;
for i = 1:length(meta_struct)
    orig_traces = meta_struct(i).raw_ORIG;
    ml_traces = meta_struct(i).raw_ML;
    orig_nc = [orig_traces.ncID];
    orig_set = [orig_traces.setID];
    for j = 1:length(ml_traces)
        nucleus = ml_traces(j).ncID;
        set = ml_traces(j).setID;
        if sum((orig_set == set).*(orig_nc == nucleus)) == 1
            t_pass = t_pass + 1; 
            orig_ind = find((orig_set == set).*(orig_nc == nucleus));
            
            of = orig_traces(orig_ind).fluo;
            ot = orig_traces(orig_ind).time;
            op = orig_traces(orig_ind).pID;
            
            mf = ml_traces(j).fluo;
            mt = ml_traces(j).time;
            mp = ml_traces(j).pID;
            
            if length(of) < 15 || length(mf) <15
                continue
            end
            if set == 7
                break
            end
            t_fig = figure('Visible','off','Position',[0 0 1024 512]);
            subplot(1,2,1)
            plot(ot, of,'-o', 'Linewidth',1.5)
            title(['Trace ' num2str(i) ' Set: ' num2str(set) ' Nuclues: ' num2str(nucleus) ...
                ' Particle: ' num2str(op) ' (Original)'])
            grid on
            axis([0 max([ot mt]) 0 1.1*max([of mf])]);
            subplot(1,2,2)
            plot(mt, mf,'-o', 'Linewidth',1.5)
            title(['Trace ' num2str(i) ' Set: ' num2str(set) ' Nuclues: ' num2str(nucleus) ...
                ' Particle: ' num2str(mp) ' (Weka)'])
            grid on
            axis([0 max([ot mt]) 0 1.1*max([of mf])]);
            saveas(t_fig, [outpath,'/traces/', 'set_' num2str(set) '_trace_' num2str(j) '.png'],'png');
        end    
    end
end
%%
spot_metric_struct = struct;
for h = 1:length(key_names)
   for k = 3:3
       p = meta_struct(k).(['particles_' key_names{h}]);
       sp = meta_struct(k).(['spots_' key_names{h}]);
       st = struct;
       for i = 1:length(p)
          %initialize the final structure for analysis
          trace_length = length([p(i).Frame]);
          st(i).Frame = zeros(1, trace_length);
          st(i).Snippet = cell(1, trace_length);
          st(i).Position = cell(1, trace_length);
%            st(i).Snippet3 = cell(1, trace_length);
          st(i).ParticleID = p(i).OriginalParticle;
          st(i).particles = p(i);
          %add in particle frame information
          st(i).Frame = [p(i).Frame];
          
          %add in particle snippets (will do a max projection, so one snippet per
           %particle per frame) and spot positions (at the brightest z-plane)
          for t = 1:trace_length
              index = p(i).Index(t);
              f = p(i).Frame(t);
              if ~isempty(sp(f).Fits)
                  fits = sp(f).Fits(index);
                  x = fits.xDoG(fits.z == fits.brightestZ);
                  y = fits.yDoG(fits.z == fits.brightestZ);
                  z = fits.brightestZ;
                  st(i).Position{t} = [x,y,z];
                  st(i).Snippet{t} = fits.Snippet(fits.z == fits.brightestZ);
              end
          end
       end
       spot_metric_struct(k).(['spot_stats_' key_names{h}]) = st;
   end
end
%%
seven_ML = spot_metric_struct(3).spot_stats_ML;
seven_ORIG = spot_metric_struct(3).spot_stats_ORIG;

particle_ORIG = 84;
particle_ML = 2;

ML_p_stats = seven_ML([seven_ML.ParticleID]==particle_ML);
ORIG_p_stats = seven_ORIG([seven_ORIG.ParticleID]==particle_ORIG);