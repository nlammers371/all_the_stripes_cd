%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\mHMM\';
sub_paths = {'weka\', 'orig\'};
key_names = {'ML', 'ORIG'};
% prefix = 'eve2_20sec_';
project = 'eve2_ml_orig_compare';
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
position = '';
date_str = '2017-07-09-';
suffix = 'eve2_20sec_19';
fName = [date_str position suffix];
%% Load Compiled Particles and Spot Mats
particles = struct;
for i = 1:length(sub_paths)
    sp = sub_paths{i};
    kn = key_names{i};
    raw_particles = load([folder_path sp fName '/' 'CompiledParticles.mat']);
    particles.(key_names{i}) = raw_particles.CompiledParticles;
    particles.(key_names{i}) = particles.(key_names{i})([particles.(key_names{i}).nc]==14);
    raw_spots = load([folder_path sp fName '/' 'Spots.mat']);
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

        if length(of) < 30 && length(mf) <30
            continue
        end
        
        t_fig = figure('Visible','off','Position',[0 0 1024 512]);
        subplot(1,2,1)
        plot(ot, of,'-o', 'Linewidth',1.5)
        title(['Trace ' num2str(j) ' Nuc.: ' num2str(nucleus) ...
            ' Part.: ' num2str(op) ' (Orig)'])
        grid on
        axis([0 max([ot mt]) 0 1.1*max([of mf])]);
        subplot(1,2,2)
        plot(mt, mf,'-o', 'Linewidth',1.5)
        title(['Trace ' num2str(j) ' Nuclues: ' num2str(nucleus) ...
            ' Particle: ' num2str(mp) ' (Weka)'], 'FontSize',7)
        grid on
        set(gca,'fontsize',10)
        axis([0 max([ot mt]) 0 1.1*max([of mf])]);
        saveas(t_fig, [outpath,'/traces/', 'set_' fName '_trace_' num2str(j) '.png'],'png');
    end    
end
% %%
% sm_struct = struct;
% for h = 1:length(key_names)
%    p = particles.(key_names{h});
%    sp = spots.(key_names{h});
%    st = struct;
%    for i = 1:length(p)
%       %initialize the final structure for analysis
%       trace_length = length([p(i).Frame]);
%       st(i).Frame = zeros(1, trace_length);
%       st(i).Snippet = cell(1, trace_length);
%       st(i).Position = cell(1, trace_length);
%       st(i).CentralIntensity = cell(1, trace_length);
%       st(i).GaussianIntensity = cell(1, trace_length);
%       
%       st(i).ParticleID = p(i).OriginalParticle;
%       st(i).particles = p(i);
%       %add in particle frame information
%       st(i).Frame = [p(i).Frame];
% 
%       %add in particle snippets (will do a max projection, so one snippet per
%        %particle per frame) and spot positions (at the brightest z-plane)
%       for t = 1:trace_length
%           index = p(i).Index(t);
%           f = p(i).Frame(t);
%           if ~isempty(sp(f).Fits)
%               fits = sp(f).Fits(index);
%               x = fits.xDoG(fits.z == fits.brightestZ);
%               y = fits.yDoG(fits.z == fits.brightestZ);
%               z = fits.brightestZ;
%               st(i).Position{t} = [x,y,z];
%               st(i).Snippet{t} = fits.Snippet(fits.z == fits.brightestZ);
%               st(i).CentralIntensity{t} = fits.CentralIntensity(fits.z == fits.brightestZ);
%               st(i).GaussianIntensity{t} = fits.GaussianIntensity(fits.z == fits.brightestZ);
%           end
%       end
%    end
%    sm_struct.(['spot_stats_' key_names{h}]) = st;
% end
% save([outpath,'/', 'set_' num2str(set_num) '_spot_info.mat'],'sm_struct')
%% AP Histograms
% vector of unique ap positions
ap_vec = unique(floor(100*[[particles.ORIG.APpos] [particles.ML.APpos]]));
ap_counts = zeros(2,length(ap_vec));
for j = 1:length(key_names)
    kn = key_names{j};
    for i = 1:length(particles.(kn))
        ap = floor(100*mean(particles.(kn)(i).APpos));
        ap_counts(j,ap_vec==ap) = ap_counts(j,ap_vec==ap) + length(particles.(kn)(i).Fluo);
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
legend(key_names{:}, 'Location', 'northeast');
grid on
axis([min(ap_vec) max(ap_vec) 0 1.2*max(max(ap_counts))]);
title(['Data Points by AP Region:' date_str position suffix]);
xlabel('AP Region (%)');
saveas(hist_fig, [outpath,'/', 'set_' date_str position suffix '_spot_info.png'],'png')
%% Trace Length Distribution

trace_fig = figure;
hold on
orig_lengths = [];
weka_lengths = [];
for i = 1:length(particles.ORIG)
    orig_lengths = [orig_lengths length(particles.ORIG(i).Fluo)];
end
for i = 1:length(particles.ML)
    weka_lengths = [weka_lengths length(particles.ML(i).Fluo)];
end
histogram(weka_lengths,30)
histogram(orig_lengths,30)
