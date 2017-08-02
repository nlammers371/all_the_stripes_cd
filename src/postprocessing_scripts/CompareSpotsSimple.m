%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\Eve7Stripes\';
sub_paths = {'weka\', 'orig\'};
key_names = {'ML', 'ORIG'};
prefix = 'eve2_20sec_';
project = 'BACeve_ml_orig_compare';
outpath = [folder_path '/figures/' project '/' ];
% Keyword to ensure only sets from current project are pulled
if exist(outpath) ~= 7
    mkdir(outpath);
end
if exist([outpath 'traces']) ~= 7
    mkdir([outpath 'traces']);
end
%-------------Designate Which Data Set to Use for Comparison--------------%
% set_num = 9;
position = '150_umP_';
date_str = '2017-05-14-';
suffix = '_eve_5uW';
fName = [data_str position suffix];
%% Load Compiled Particles and Spot Mats
particles = struct;
for i = 1:length(sub_paths)
    sp = sub_paths{i};
    kn = key_names{i};
    raw_particles = load([folder_path sp fName '/' 'CompiledParticles.mat']);
    particles.(key_names{i}) = raw_particles.CompiledParticles;
    raw_spots = load([folder_path sp date_str prefix num2str(set_num) '/' 'Spots.mat']);
    spots.(key_names{i}) = raw_spots.Spots;
end
%%
t_pass = 1;
O_particles = particles.ORIG;
orig_nc = [O_particles.Nucleus];
% ml_traces = particles.ML.fluo;
for j = 1:length(particles.ML)
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
        saveas(t_fig, [outpath,'/traces/', 'set_' num2str(set_num) '_trace_' num2str(j) '.png'],'png');
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
% vector of unique ap positions
ap_vec = unique(floor([[particles.ORIG.APpos] [particles.ML.APpos]]));
ap_counts = zeros(2,length(ap_vec));
for j = 1:length(key_names)
    kn = key_names{j};
    for i = 1:length(particles.(kn))
        ap = floor(mean(particles.(kn)(i).APpos));
        ap_counts(i,ap_vec==ap) = ap_counts(1,ap_vec==ap) + length(particles.(kn)(i).fluo);
    end
end

hist_fig = figure;
hold on

colormap('jet');
cm = colormap;

bar(ap_vec, ap_counts(1,:), 'Facecolor',cm(5,:),...
    'EdgeColor','black','FaceAlpha',.5,'BarWidth', 1);
bar(ap_vec, ap_counts(2,:), 'Facecolor',cm(25,:),...
    'EdgeColor','black','FaceAlpha',.5,'BarWidth', 1);
legend('Original Pipeline', 'Weka');

title('Data Points by AP Region');
xlabel('AP Region (%)');
saveas(hist_fig, [outpath,'/', 'set_' num2str(set_num) '_spot_info.png'],'png')