addpath('../utilities/');

date_str = '2017_06_21';
meta_folder = 'inference_results';
folder_path = ['../../inference_results/' date_str '/'];
files = dir(folder_path);

filenames = {};
for i = 1:length(files)
    if strfind(files(i).name,'.mat') > 0
        filenames = [filenames {files(i).name}];
    end
end

outpath = ['../../inference_results/' date_str '/summaries/'];
if (exist(outpath, 'dir') ~= 7)
    mkdir(outpath);
end
%Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'
    load(['../../inference_results/' date_str '/' filenames{f}], 'output');
    for fn = fieldnames(output)'
        glb_all(f).(fn{1}) = output.(fn{1});
    end
end

OutPath = ['../../fig/' 'experimental_system' '/'];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

%memory assumed for inference
w = 8;

%states used for final inference
K = 3;

%time step used in experimental data
dT = 18.75;

%!!!set AP range
ap_range = [31:54];

glb_all = glb_all(ismember(floor([glb_all.AP]),ap_range));

%load traces (saved as "preprocessed")
traces_all = load(['../../processed_data/eveSet_2017_06_15.mat']);
traces_all = traces_all.interp_struct;
traces_all = traces_all(ismember(floor([traces_all.AP]),ap_range));


%% Generate Occupancy-Loading Rate Plot
%filter for desired AP range
glb_ap_index = [glb_all.AP];
alpha = glb_all(1).alpha;
occupancy_ap = zeros(K,length(glb_all));
occupancy_ap_v = occupancy_ap;
emission_rates = zeros(K,length(glb_all));
for i = 1:length(glb_all)
    [emission_rates(:,i), ranked_v] = sort([glb_all(i).v]);
    A_log = reshape(glb_all(i).A_log,K,K);
    A_log = A_log(ranked_v, ranked_v);
    noise = glb_all(i).noise;
    pi0_log = glb_all(i).pi0_log;
    pi0_log = pi0_log(ranked_v);
                                   
    [V,D] = eig(exp(A_log));
    steady = V(:,1)./sum(V(:,1));
    occupancy_ap(:,i) = steady;
end

%create bubble chart conveying dwell times and emission rates by AP
bubble_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(50,:) ; cmap(30,:); cmap(20,:)]);

init_vec = reshape(emission_rates,1,[]);
ap_vec = repelem([glb_all.AP],1,K);
c_mat = repmat(cm,length([glb_all.AP]),1);
hold on

plot(reshape(ap_vec,K,[])', emission_rates','black','LineWidth',1);

%h = colorbar;

% shapes = {'o','^','s'};

sctr_emission = scatter(ap_vec', init_vec , 75, c_mat, 'o','filled', 'MarkerEdgeColor', 'black') ;

% h;

box off
axis([(min(ap_range)-1) (max(ap_range)+1) -500 1.2*max(init_vec)])

set(gca,'FontName','Lucida Sans Regular')
grid on
title('State Fluorescence by AP Position');
xlabel('AP Position (%)');
ylabel('Fluorescence (A.U)');
% ylabel(h, 'Average Occupancy')

saveas(bubble_fig, [OutPath '/emission.eps'], 'epsc');
saveas(bubble_fig, [OutPath '/emission.png'], 'png');
hold off

%Occupancy
occ_fig = figure;
hold on
plot(reshape(ap_vec,K,[])', occupancy_ap','black','LineWidth',1);

sctr_occ = scatter(ap_vec', reshape(occupancy_ap,1,[])' , 75, c_mat, 'o','filled', 'MarkerEdgeColor', 'black') ;

% h;

box off
axis([(min(ap_range)-1) (max(ap_range)+1) 0 1.2*max(max(occupancy_ap))])

set(gca,'FontName','Lucida Sans Regular')

title('State Occupancy by AP Position');
ylabel('Occupancy Share');
xlabel('AP Position (%)');

grid on

saveas(occ_fig, [ OutPath '/occupancy.eps'], 'epsc');
saveas(occ_fig, [ OutPath '/occupancy.png'], 'png');

%% Dwell Times
dwell_times = zeros(K,length(glb_all));

for i = 1:length(glb_all)
    Rcol = glb_all(i).R;
%     R = reshape(Rcol,K,K);
    [~, ranked_v] = sort(glb_all(i).v);
    A = reshape(glb_all(i).A,K,K);
    A = A(ranked_v, ranked_v);
    R = reshape(glb_all(i).R,K,K);
    glb_all(i).R = reshape(R(ranked_v, ranked_v),1,[]);
    if ~isreal(Rcol)||(sum(Rcol<0)>K)
        out = prob_to_rate_fit_sym(A, dT, 'gen', .005, 1);            
        glb_all(i).R_fit = out.R_out;
        r_diag = diag(out.R_out);
    else
        glb_all(i).R_fit = glb_all(i).R;
        r_diag = diag(reshape(glb_all(i).R,K,K)); 
    end
    dwell_times(:,i) = -r_diag.^-1;
end
%%
dwell_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(50,:) ; cmap(30,:); cmap(20,:)]);

ap_vec = repelem([glb_all.AP],1,K);

%!!! Need to update to incorporate real data
hold on
error_vec = 2 + 10*rand(1, length(ap_vec));
%errorbar(ap_vec,reshape(dwell_times,1,[]),error_vec,'.k')
c_mat = repmat(cm,length([glb_all.AP]),1);

% shapes = {'o','^','s'};
% h = zeros(3, 1);
% for i = 1:K
%     h(i) = plot(0,0,shapes{i}, 'Color', 'black','MarkerFaceColor', cm );
% end
plot([glb_all.AP], dwell_times,'black','LineWidth',1);

sctr = scatter(ap_vec', reshape(dwell_times,1,[])', 75, c_mat, 'o', 'filled', 'MarkerEdgeColor', 'black') ;


% legend(h, 'OFF','ON1','ON2');

box off
set(gca,'FontName','Lucida Sans Regular')
axis([(min([glb_all.AP])-1) (max([glb_all.AP])+1) 0 15*ceil(max(max(dwell_times))/15)])
title('Dwell Times by AP Position');
xlabel('AP Position (%)');
ylabel('Dwell Time (s)');
grid on 
saveas(dwell_fig, [OutPath '/dwell_times.eps'], 'epsc');
saveas(dwell_fig, [OutPath '/dwell_times.png'], 'png');

%% Look at Implied k_off and k_on, burst size and burst freq
%For now I am assuming simple 2 promoter model...unclear if this is really
%justfied given strange results for 3 state inference
k_on_vec = zeros(1,length(glb_all));
k_off_vec = zeros(1,length(glb_all));
for i = 1:length(glb_all)
    R = reshape(glb_all(i).R_fit,K,K);
    % R21 = 2*k_on
    k_on_vec(i) = R(2,1) / 2;
    % R12 = k_off
    k_off_vec(i) = R(1,2);
end

k_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(35,:); cmap(10,:)]);

ap_vec = repmat([glb_all.AP],1,2);

%!!! Need to update to incorporate real data
hold on
c_mat = repelem(cm,length([glb_all.AP]),1);

plot([glb_all.AP], k_on_vec,'black','LineWidth',1);
plot([glb_all.AP], k_off_vec,'black','LineWidth',1);
% plot([0],[0],[0],[0], 'Marker', 'o', 'Color',cm(1,:));
sctr = scatter(ap_vec, [k_on_vec,k_off_vec], 75, c_mat, 'o', 'filled', 'MarkerEdgeColor', 'black') ;
% legend('k_{off}','k_{on}');
title('k_{on} (blue) and k_{off} (green) by AP Position');
xlabel('AP Position (%)');
ylabel('Rate (s^-1)');
grid on 
saveas(k_fig, [OutPath '/switching_rates.eps'], 'epsc');
saveas(k_fig, [OutPath '/switching_rates.png'], 'png');

cycle_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(45,:)]);

ap_vec = [glb_all.AP];

off_times = dwell_times(1,:);
off_occ = occupancy_ap(1,:);
t_cycle = off_times./off_occ;
%!!! Need to update to incorporate real data
hold on

plot([glb_all.AP], t_cycle.^-1,'black','LineWidth',1);

sctr = scatter(ap_vec, t_cycle.^-1, 75, cm, 'o', 'filled', 'MarkerEdgeColor', 'black') ;
% legend('k_{off}','k_{on}');
title('Characteristic Cycle Frequency by AP Position');
xlabel('AP Position (%)');
ylabel('Cycle Frequency (s^-1)');
grid on 
saveas(cycle_fig, [OutPath '/cycle_times.eps'], 'epsc');
saveas(cycle_fig, [OutPath '/cycle_times.png'], 'png');

n_fig = figure;

colormap((colormap('jet')));
cmap = colormap;
cm = flipud([cmap(5,:)]);

ap_vec = [glb_all.AP];

dp_count = [glb_all.N];
%!!! Need to update to incorporate real data
hold on

plot([glb_all.AP], dp_count,'black','LineWidth',1);

sctr = scatter(ap_vec, dp_count, 75, cm, 'o', 'filled', 'MarkerEdgeColor', 'black') ;
% legend('k_{off}','k_{on}');
title('Number of Data Points by AP Position');
xlabel('AP Position (%)');
ylabel('N Data Points');
grid on 
saveas(n_fig, [OutPath '/n_dp.eps'], 'epsc');
saveas(n_fig, [OutPath '/n_dp.png'], 'png');

%% Trace interp check
interp_fig = figure;
plot(traces_all(i).time_orig / 60, traces_all(i).fluo_orig, traces_all(i).time / 60,traces_all(i).fluo,'LineWidth', 1)
legend('Original', 'Interpolated');
title('Comparing Interpolated Trace with Original')
xlabel('Time (min)');
ylabel('AU');
saveas(interp_fig, [OutPath '/interp.eps'], 'epsc');
saveas(interp_fig, [OutPath '/interp.png'], 'png');


%% Viterbi Plots
%states
K = 3;
%memory
w = 8;
%alpha (lenght of MS2 Loops in time steps)
alpha = glb_all(1).alpha;

%subgroup to summarize
gr = 40.5;

colormap('parula');
cmap = colormap;
cm1 = cmap(20,:);
% warning(ms2_loading_coeff (alpha, w))
%viterbi plots
for i = 9
    tic
    fluo = traces_all(i).fluo;
    gb_fit = glb_all((1*(ceil([glb_all.AP]) == traces_all(i).AP) + 1*(floor([glb_all.AP]) == traces_all(i).AP))==1 );
    v_fit = viterbi (fluo, gb_fit.v, gb_fit.noise, gb_fit.pi0_log, ...
                                     gb_fit.A_log, K, w, alpha);
    toc
    v_fluo = v_fit.fluo_viterbi;
end
%%
%     naivePadded = [zeros(1,w-1) v(naive)];
%     register = zeros(1,w);
%     wt_vec = ones(1,w);
%     for m = 1:w
%         if m < alpha 
%             factor = .5*(m-1)/alpha + .5*m/alpha;
%             wt_vec(m) = factor;
%         elseif m > alpha && m < alpha + 1
%             factor = (alpha-m+1)*.5*(1+(m-1)/alpha) + m-alpha;
%             wt_vec(m) = factor;
%         end
%         
%     end
%     wt_vec = fliplr(wt_vec);
    F = zeros(1,length(naive));
    for i = w:(w-1+length(naive))
        register = naivePadded(i-w+1:i);
        F(i-w+1) = sum(register.*wt_vec);
    end       
       
	% Create a filename.
	traceFName = sprintf('trace_%d.eps', i);
    stateFName = sprintf('state_%d.eps', i);
	traceFileName = fullfile(OutPath, traceFName);
	stateFileName = fullfile(OutPath, stateFName);	
       
    % plot true vs viterbi fluorescence values
    figure('Position', [100, 100, 500, 300],'Visible', 'on');
    hold on
    f = plot(fluo, 'LineWidth', 3, 'Color', 'k');
    plot(F, 'LineWidth', 3, 'Color', cm1);
    ylim([0, 25000]);
    ax = ancestor(f, 'axes');
    xrule = ax.XAxis;
    xrule.FontSize = 18;
    yrule = ax.YAxis;
    yrule.FontSize = 18;
    % Then save it
	saveas(f, traceFileName, 'epsc');
    hold off
    
    %naive state
    figure('Position', [100, 100, 500, 300],'Visible', 'off');
    hold on
    f = stairs(v(naive), 'black', 'LineWidth', 2);
    ax = ancestor(f, 'axes');
    title(['Example Trace with Viterbi Fit (' num2str(gr) ')']);
    xrule = ax.XAxis;
    xrule.FontSize = 18;
    yrule = ax.YAxis;
    yrule.FontSize = 18;
    set(gca,'FontName','Lucida Sans Regular')
%     yticks([0, 20, 40, 60, 80, 100]);
    saveas(f, stateFileName, 'eps');
    hold off
end