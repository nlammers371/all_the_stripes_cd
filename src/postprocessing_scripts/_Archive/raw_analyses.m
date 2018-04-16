%Script to run analyses on inference data
addpath('../utilities/');

%------------------------------Set System Params--------------------------%
%memory assumed for inference
w = 7;
%Tres
dT = 20;

datatype = 'weka';
date_str = '2017-09-25';
project = ['eve7stripes_inf_2017_09_25'];

OutPath = ['../../fig/experimental_system' '/' project '/' date_str '/other_analyses/' ];
if exist(OutPath) ~= 7
    mkdir(OutPath);
end

%load traces (saved as "interp_struct")
traces_all = load(['../../dat/' project '/inference_traces_t' num2str(dT) '_' project '.mat']);
traces_all = traces_all.interp_struct;

%% Noise vs Mean Scatter
set_list = unique([traces_all.setID]);
%set id vec
set_vec = [];
%stripe id vec
stripe_vec = [];
%mean fluo vec
fluo_vec = [];
%noise vec 
noise_vec = [];
%loop through the sets
for s = 1:length(set_list)
    set_curr = set_list(s);
    set_traces = traces_all([traces_all.setID]==set_curr);
    %Fliter for only traces in stripe centers
    set_traces = set_traces([set_traces.stripe_sub_id]==0);
    stripes = unique([set_traces.stripe_id]);
    for st = stripes
        st_set = set_traces([set_traces.stripe_id]==st);
        %store info about cumulative fluroescence
        cf_vec = [];
        for i = 1:length(st_set)
            cf_vec = [cf_vec sum(st_set(i).fluo)];
        end
        fluo_vec = [fluo_vec mean(cf_vec)];
        noise_vec = [noise_vec std(cf_vec)];
        stripe_vec = [stripe_vec st];
        set_vec = [set_vec set_curr];
    end
end
sym_vec = {'s','o','d','p','h','*','+'};
fano_fig = figure;
colormap('jet');
cm = colormap;
inc = floor(size(cm,1)/length(set_list));
cmSet = cm(1:length(set_list)*inc,:);
set_ind = [];
for i = set_vec
    set_ind = [set_ind find(set_list==i,1)];
end
hold on
for st = 1:7
    psym = sym_vec{st};
    si = set_ind(stripe_vec==st);
    a = scatter(fluo_vec(stripe_vec==st),noise_vec(stripe_vec==st),75,...
        cm(inc*si,:),psym);
    set(a,'LineWidth',2);
end
grid on
title('Mean Cumulative Fluorescence vs Noise');
xlabel('Mean Cumulative (AU)')
ylabel('Standard Deviation (AU)')

saveas(fano_fig, [OutPath '/fano_fig.eps'], 'epsc');
saveas(fano_fig, [OutPath '/fano_fig.png'], 'png');
%%
stripe_id = 1:7;
md = [];
vd = [];
n_boots = 100;
for s = 1:7
    f_diffs = [];
    s_traces = traces_all([traces_all.stripe_id]==s);
    for i = 1:length(s_traces)
        f_diffs = [f_diffs abs(diff(s_traces(i).fluo))];
    end
    md_sub = [];
    for b = 1:n_boots
        bs = randsample(f_diffs,length(f_diffs),true);
        md_sub = [md_sub mean(bs)];
    end
    md = [md mean(md_sub)];
    vd = [vd std(md_sub)];
%     histogram(f_diffs)
end

% grid on
jump_fig = figure;
colormap('winter')
cm = colormap;
hold on
plot(1:7,md,'LineWidth',2,'Color',cm(30,:))
scatter(1:7,md,75,'MarkerFaceColor',cm(30,:),'MarkerEdgeColor','black','LineWidth',1.5)
grid on
xlabel('Eve Stripe')
ylabel('Average Jump Size')
title('Average Jump Size by Stripe')

saveas(jump_fig, [OutPath '/jump_fig.eps'], 'epsc');
saveas(jump_fig, [OutPath '/jump_fig.png'], 'png');

%% -----------------------Line Scan Test-------------------------------- %%
% set granularity
granularity = 25000;
% generate structures to store results
max_fluo = max([traces_all.fluo]);
n_bins = ceil(max_fluo/granularity);
slice_vec = granularity:granularity:n_bins*granularity;
slice_struct = struct;
for s = 1:7
    slice_struct(s).gap_cell = cell(1,n_bins);
end
% Look for cross points
for tr = 1:length(traces_all)
    fluo = traces_all(tr).fluo;
    time = traces_all(tr).time;
    stripe = traces_all(tr).stripe_id;
    gap_cell = slice_struct(stripe).gap_cell;
    for cr = 1:n_bins
        cross = slice_vec(cr);
        fluo_bin = fluo > cross;
        dfb = diff(fluo_bin);
        gaps = diff(time([0 dfb]~=0));
        gaps = gaps(gaps>100);
        gc = gap_cell{cr};
        gap_cell{cr} = [gc gaps];        
    end
    slice_struct(stripe).gap_cell = gap_cell;
end
%% 
gc = slice_struct(4).gap_cell;
gap_fig = figure;
hold on
for n = 1:n_bins
    gaps = gc{n};
    plot(linspace(slice_vec(n),slice_vec(n),length(gaps)), gaps,'o')
    plot(slice_vec(n),mean(gaps),'s','MarkerSize',7,'LineWidth',2)
end
    
