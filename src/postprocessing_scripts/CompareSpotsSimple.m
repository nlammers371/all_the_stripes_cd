%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\DropboxSingleTraces\';
sub_paths = {'Eve2_ML/', 'Eve2_orig/'};
key_names = {'ML', 'ORIG'};
prefix = 'eve2_20sec_';
project = 'mHMMeve2_ml_comparisons';
outpath = ['../projects/' project '/' ];
% Keyword to ensure only sets from current project are pulled
if exist(outpath) ~= 7
    mkdir(outpath);
end
%% Load Compiled Particles and Spot Mats
particles = struct;
set_num = 7;
date_str = '2017-05-30-';
for i = 1:length(sub_paths)
    sp = sub_paths{i};
    kn = key_names{i};
    raw_particles = load([folder_path sp date_str prefix num2str(set_num) '/' 'CompiledParticles.mat']);
    particles.(key_names{i}) = raw_particles.CompiledParticles;
    raw_spots = load([folder_path sp date_str prefix num2str(set_num) '/' 'Spots.mat']);
    spots.(key_names{i}) = raw_spots.Spots;
end
%%
t_pass = 1;
O_particles = particles.ORIG;
orig_nc = [O_particles.Nucleus];
for j = 1:length(ml_traces)
    nucleus = particles.ML(j).Nucleus;
    if sum(orig_nc == nucleus) == 1
        t_pass = t_pass + 1; 
        orig_ind = find(orig_nc == nucleus);

        of = particles.ORIG(orig_ind).Fluo;
        ot = particles.ORIG(orig_ind).Frame;
        op = particles.ORIG(orig_ind).OriginalParticle;

        mf = particles.ML(j).Fluo;
        mt = particles.ML(j).Frame;
        mp = particles.ML(j).OriginalParticle;

        if length(of) < 15 || length(mf) <15
            continue
        end
        if set == 7
            break
        end
        t_fig = figure('Visible','off','Position',[0 0 1024 512]);
        subplot(1,2,1)
        plot(ot, of,'-o', 'Linewidth',1.5)
        title(['Trace ' num2str(j) ' Set: ' num2str(set_num) ' Nuclues: ' num2str(nucleus) ...
            ' Particle: ' num2str(op) ' (Original)'])
        grid on
        axis([0 max([ot mt]) 0 1.1*max([of mf])]);
        subplot(1,2,2)
        plot(mt, mf,'-o', 'Linewidth',1.5)
        title(['Trace ' num2str(j) ' Set: ' num2str(set_num) ' Nuclues: ' num2str(nucleus) ...
            ' Particle: ' num2str(mp) ' (Weka)'])
        grid on
        axis([0 max([ot mt]) 0 1.1*max([of mf])]);
        saveas(t_fig, [outpath,'/traces_new/', 'set_' num2str(set_num) '_trace_' num2str(j) '.png'],'png');
    end    
end
%%
sm_struct = struct;
for h = 1:length(key_names)
   p = particles.(key_names{h});
   sp = spots.(key_names{h});
   st = struct;
   for i = 1:length(p)
      %initialize the final structure for analysis
      trace_length = length([p(i).Frame]);
      st(i).Frame = zeros(1, trace_length);
      st(i).Snippet = cell(1, trace_length);
      st(i).Position = cell(1, trace_length);
      st(i).CentralIntensity = cell(1, trace_length);
      st(i).GaussianIntensity = cell(1, trace_length);
      
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
              st(i).CentralIntensity{t} = fits.CentralIntensity(fits.z == fits.brightestZ);
              st(i).GaussianIntensity{t} = fits.GaussianIntensity(fits.z == fits.brightestZ);
          end
      end
   end
   sm_struct.(['spot_stats_' key_names{h}]) = st;
end
save([outpath,'/', 'set_' num2str(set_num) '_spot_info.mat'],'sm_struct')
%%
seven_ML = sm_struct(1).spot_stats_ML;
seven_ORIG = sm_struct(1).spot_stats_ORIG;

particle_ORIG = 63;
particle_ML = 63;

ML_p_stats = seven_ML([sm_struct.spot_stats_ML.ParticleID]==particle_ML);
ORIG_p_stats = seven_ORIG([seven_ORIG.ParticleID]==particle_ORIG);
