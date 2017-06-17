load('C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\data\CompiledParticles_2017_02_08_eve23_4.mat')

%Filter for NC14 Traces
cp_14 = CompiledParticles([CompiledParticles.nc]==14);

%% Spatial Distribution of Traces
histogram([cp_14.APPos])
% Looks like we have two stripes represented here

%% Examine Raw Traces
for i = 1:length(cp_14)
    hold on
    plot(cp_14(i).Fluo)
    title(['Trace ' num2str(i)])
    
end
hold off
%Some remarkably well-distinguished Peaks in there
%% Look at Distribution of Deltas 
delta_rates = [];

for i = 1:length(cp_14)
    diff_rates = diff(cp_14(i).Fluo);
    delta_rates = [delta_rates diff_rates];
end
histogram(abs(delta_rates));

%% Fluo Distribution
histogram([cp_14.Fluo]);
%% Bivariate Hist
delta_rates = [];
f_levels = [];
for i = 1:length(cp_14)
    f = cp_14(i).Fluo;
    diff_rates = diff(f) / 6;
    delta_rates = [delta_rates diff_rates];
    mean_fluo = (f(1:end-1) + f(2:end)) / 2;
    f_levels = [f_levels mean_fluo];
end
histogram2(f_levels', delta_rates', 'FaceColor','flat');
%% Autocorrelation...
%Filter for nc14 trecs
fluo_vec_14 = AllTracesVector(:,ncFilter(:,4));
fluo_vec_14(isnan(fluo_vec_14)) = 0;
n_lag = 50;
tr_ct = 0;
gaussFilter = gausswin(5);
gaussFilter = gaussFilter/sum(gaussFilter);
for i = 1:size(fluo_vec_14,2)
    f_vec = fluo_vec_14(:,i);
    clipped_fluo = f_vec(find(f_vec,1):find(f_vec,1,'last'));
    if length(clipped_fluo) >= n_lag + 1
        tr_smooth = conv(clipped_fluo,gaussFilter,'same');
        acf_smooth = autocorr(tr_smooth, n_lag);
        acf = autocorr(clipped_fluo, n_lag);
        if tr_ct == 0
            acf_mat_smooth = acf_smooth';
            acf_mat_dd_smooth = diff(acf_smooth,2)';
            acf_mat = acf';
            acf_mat_dd = diff(acf,2)';
        else
            acf_mat_smooth = vertcat(acf_mat_smooth,acf_smooth');
            acf_mat_dd_smooth = vertcat(acf_mat_dd_smooth,diff(acf_smooth,2)');
            acf_mat = vertcat(acf_mat,acf');
            acf_mat_dd = vertcat(acf_mat_dd,diff(acf,2)');
        end
        tr_ct = tr_ct + 1;
    end
end
colormap('winter');
cm = colormap;
mean_acf = mean(acf_mat);
mean_acf_dd = mean(acf_mat_dd);
mean_acf_smooth = mean(acf_mat_smooth);
mean_acf_dd_smooth = mean(acf_mat_dd_smooth);
% std_acf = std(acf_mat);
autoFig = figure;
hold on
plot(0:(n_lag),mean_acf,0:(n_lag),mean_acf_smooth, 'LineStyle', '-','Marker', 'o', 'MarkerSize', 2, 'LineWidth',1.5)
legend('Raw Traces', 'Smoothed Traces');
% errorbar(0:6:(6*n_lag),mean_acf, std_acf, 'Color','black')
% plot(linspace(0,6*n_lag),linspace(0,0),'LineStyle','--','Color',cm(40,:))
xlabel('Lag (s)');
ylabel('Autocorrelation');
title('Mean Autocorrelation Function For NC 14 Traces');
grid on
% axis([0 6*n_lag -.4 1.0])
hold off

autoFigDD = figure;
hold on
plot(2:(n_lag),mean_acf_dd,2:(n_lag),mean_acf_dd_smooth, 'LineStyle', '-','Marker', 'o', 'MarkerSize', 2, 'LineWidth',1.5)
legend('Raw Traces', 'Smoothed Traces');
% errorbar(0:6:(6*n_lag),mean_acf, std_acf, 'Color','black')
% plot(linspace(0,6*n_lag),linspace(0,0),'LineStyle','--','Color',cm(40,:))
xlabel('Lag (time steps)');
title('Mean Second Derivative of Autocorrelation For NC 14 Traces');
plot(linspace(0,n_lag+1),linspace(0,0),'LineStyle','--','Color','black')
axis([0 n_lag -.015 .02])
grid on
hold off

%% Make Histograms
dd_max_vec = zeros(1,length(acf_mat_dd));
dd_max_vec_smooth = zeros(1,length(acf_mat_dd));
for i = 1:length(acf_mat_dd)
    [m, m_ind] = max([zeros(1,6), acf_mat_dd(i,6:end)]);
    dd_max_vec(i) = m_ind;
    [ms, ms_ind] = max([zeros(1,6), acf_mat_dd_smooth(i,6:end)]);
    dd_max_vec_smooth(i) = ms_ind;
end
histFig = figure;
hold on
histogram(dd_max_vec)
histogram(dd_max_vec_smooth);
legend('Raw Traces', 'Smoothed Traces');
title('Histogram of Maxima: Second Derivative of Trace Autocorrelation');

hold off

%% Estimate Noise

%Filter for nc14 trecs
fluo_vec_14 = AllTracesVector(:,ncFilter(:,4));
fluo_vec_14(isnan(fluo_vec_14)) = 0;
abs_diffs = [];
max_fluo = 0;
for i = 1:size(fluo_vec_14,2)
    f_vec = fluo_vec_14(:,i);
    abs_diff = abs(diff(f_vec(find(f_vec,1):find(f_vec,1,'last'))));
    abs_diffs = [abs_diffs abs_diff'];
    max_fluo = max(max_fluo,max(f_vec));
end
display(['Mean diff: ' num2str(mean(abs_diffs))])
display(['Sig-to-Noise ' num2str(max_fluo/mean(abs_diffs))])    

