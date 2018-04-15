addpath('../utilities');
%-------------------------- Analyze Viterbi Fits -------------------------%

date_str = '2017-09-25';
project = 'eve7stripes_inf_2017_09_25';
%Memory
w = 7;
K = 3;
%Use aggregate parameter fits?
aggregate_fits = 0;
%If one, single trace plots will be generated
trace_plots = 0;
%State count
if aggregate_fits 
    ViterbiPath = ['../../out/' project '/' date_str '/mem' ...
        num2str(w) '_states' num2str(K) '/viterbi_fits_aggregate/' ];
    out_dir = ['..\..\fig\experimental_system\' project '\' date_str '\mem' ...
    num2str(w) '_states' num2str(K) '_viterbi_analysis_aggregate\'];
else    
    ViterbiPath = ['../../out/' project '/' date_str '/mem' ...
        num2str(w) '_states' num2str(K) '/viterbi_fits/' ];
    out_dir = ['..\..\fig\experimental_system\' project '\' date_str '\mem' ...
    num2str(w) '_states' num2str(K) '_viterbi_analysis\'];
end


ViterbiFiles = dir(ViterbiPath);

%If necessary, make output directory
if exist(out_dir) ~= 7
    mkdir(out_dir);
end
trace_dir = [out_dir '/traces/'];
if exist(trace_dir) ~= 7
    mkdir(trace_dir);
end
t_start = 5;
%list of AP values to check
stripe_range = [1:7];
%other inference parameters
granularity = 10;
dT = 20;
max_time = 60;
b_ct = 1;
%Three State Results
v_fit_struct = struct;
%store rates and transition probs to compare across stripes
stripe_rate_array = zeros(K,K,length(stripe_range));

for bin = stripe_range
    %Track Parameters
    occupancy_array = zeros(K,max_time);
    dwell_times = cell(1,K);
    switch_counts = zeros(K,K,max_time) + 10e-6;
    F_counts = cell(1,max_time);
    Fluo_values = cell(1,max_time);
    for i = 1:length(ViterbiFiles)              
        if ~isempty(strfind(ViterbiFiles(i).name,['bin' num2str(bin)]))
            load([ViterbiPath ViterbiFiles(i).name]);
            %Extract fit parameters
            naive_path = v_fit.z_viterbi;            
            v = v_fit.fluo_viterbi;
            f_exp = v_fit.fluo_exp;
            t_exp = v_fit.time_exp;            
            for m = 1:length(t_exp)
                minute = ceil(t_exp(m)/60 + .01);                                                
                s3 = naive_path(m);                
                f_ct = zeros(1,K);
                f_val = f_exp(m);
                F_t = Fluo_values{minute};
                Fluo_values{minute} = [F_t f_val];
                s_register = naive_path(max(1,m-w+1):m);
                s_register = [s_register ones(1,w-m-1)];
                for k = 1:K
                    f_ct(k) = sum(s_register == k);
                end
                F_ct = F_counts{minute};
                F_counts{minute} = [F_ct ; f_ct];                
                %Register state occupancy
                occupancy_array(s3,minute) = occupancy_array(s3,minute) + 1;
                %Register switch counts
                if m > 1                    
                    switch_counts(s_old,s3,minute) = switch_counts(s_old,s3,minute) + 1;
                end
                s_old = s3;
            end
            if trace_plots
                cp_fig = figure('Position', [0 0 1024 1024], 'Visible', 'off');            
                subplot(2,1,1);
                hold on
                plot(f_exp, 'Linewidth', 1.5);
                plot(v, 'Linewidth', 1.5);
                legend('Raw Trace', 'Viterbi Fit');
                ylabel('AU')
                xlabel('Time Steps');
                title(['Viterbi Fit for Experimental Stripe ' num2str(bin) ' Trace ' num2str(i) ' (' num2str(K) ' States)']);
                grid on
                subplot(2,1,2);
                stairs(naive_path,'Linewidth',1.5)
                grid on
                title(['Inferred Naive Trajectory'])
                saveas(cp_fig, [trace_dir 'bin_' num2str(bin) '_trace_' num2str(i) '.png'], 'png')            
            end
            %Add dwell times for each state. Ignore trace boundaries             
            for k = 1:length(dwell_times)
                d_times = dwell_times{k};
                binary = naive_path == k;
                d_binary = diff(binary);
                starts = find(d_binary == 1) + 1;
                stops = find(d_binary == -1);
                stops = stops(stops>=min(starts));
                starts = starts(starts<=max(stops));
                gaps = stops - starts + 1;
                dwell_times{k} = [d_times gaps];
            end
        end
    end
    close all
    wt_three = zeros(1,K);
    wait_fig = figure;
    hold on        
    for k = 1:K
        subplot(K,1,k);
        hold on
        wait_bins = 1:ceil(10*60/dT);
        wt_ct = histc(dwell_times{k}, wait_bins);        
        wt_mean = sum(wait_bins.*wt_ct)/sum(wt_ct);
        wt_three(k) = wt_mean;
        wt_ct = wt_ct / sum(wt_ct);
        fit = exppdf(wait_bins,wt_mean);
        bar(wait_bins*dT, wt_ct,'BarWidth', 1,'FaceAlpha',0.5);
        plot(wait_bins*dT,fit,'LineWidth',2)        
        title(strvcat(['Distribution of Dwell Times with Poisson Fit: State ' num2str(k-1)],...
                    ['']))             
        grid on        
    end
    saveas(wait_fig, [out_dir 'bin_' num2str(bin) '_' num2str(K) 'state_fit.png'], 'png')       
    
    %Calculate Parameters of interest for region
    occupancy_array = occupancy_array ./ repmat(sum(occupancy_array),K,1);    
    %calculate R(t)
    A_array = NaN(K,K,max_time);
    R_array = NaN(size(A_array));
    window = floor(granularity/2);
    for t = t_start:max_time
        A = sum(switch_counts(:,:,max(t-window,1):min(t+window,max_time)),3);
        A = A ./ repmat(sum(A),K,1);
        A_array(:,:,t) = A;
        R = prob_to_rate(A,dT);
        if ~isreal(R)||(sum(R(:)<0)>K)
            out = prob_to_rate_fit_sym(A, dT, 'gen', .005, 1);            
            R_fit = out.R_out;
            R_array(:,:,t) = R_fit;
        else
            R_array(:,:,t) = R;
        end
    end
    %calculate v(t)
    v_array = zeros(size(occupancy_array));
    for t = t_start:max_time        
        fluo = [Fluo_values{max(t-window,1):min(max_time,t+window)}];
        F_ct = vertcat(F_counts{max(t-window,1):min(max_time,t+window)});
        F_square = zeros(K,K);
        b_agg = zeros(1,K);
        if ~isempty(F_ct)
            for k = 1:K
                F_square(k,:) = sum(F_ct.*repmat(F_ct(:,k),1,K));
                b_agg(k) = sum(fluo'.*F_ct(:,k));
            end
            % solve system                
            v_array(:,t) = linsolve(F_square,b_agg')'; 
        end
    end
    close all    
    
    v_fig = figure('Position', [0 0 1024 512]);           
    hold on
    legend_string = [];
    for k = 1:K
        plot(1:max_time, v_array(k,:), 'LineWidth',2)
        legend_string = [legend_string {['state ' num2str(k)]}];
    end
    legend(legend_string{:})
    title(['Emission Rates Shares Over Time Stripe ' num2str(bin)])        
    grid on
    saveas(v_fig, [out_dir 'stripe_' num2str(bin) '_' num2str(K) 'state_v_trends.png'], 'png')       
    
    A_fig = figure('Position', [0 0 1024 512]); 
    hold on
    legend_string = [];
    for k = 1:K^2
        row = ceil(k/K);
        col = k - K*(row-1);
        if row == col
            continue
        end
        A_vec = reshape(A_array(row,col,t_start:50),1,[]);
        plot(t_start:50, A_vec, 'LineWidth',2)
        legend_string = [legend_string {['A' num2str(row) num2str(col)]}];
    end    
    legend(legend_string{:})
    title(['Transition Probabilities Over Time Stripe ' num2str(bin)])        
    grid on
    saveas(A_fig, [out_dir 'stripe_' num2str(bin) '_' num2str(K) 'state_A_trends.png'], 'png')       
    
    R_fig = figure('Position', [0 0 1024 512]); 
    hold on
    legend_string = [];
    for k = 1:K^2
        row = ceil(k/K);
        col = k - K*(row-1);
        if row == col
            continue
        end
        R_vec = reshape(R_array(row,col,t_start:50),1,[]);
        plot(t_start:50, R_vec, 'LineWidth',2)
        legend_string = [legend_string {['R' num2str(row) num2str(col)]}];
    end
    legend(legend_string{:})
    title(['Transition Rates Over Time Stripe ' num2str(bin)])        
    grid on
    saveas(R_fig, [out_dir 'stripe_' num2str(bin) '_' num2str(K) 'state_R_trends.png'], 'png')       
    
    
    occupancy_fig = figure('Position', [0 0 1024 512]);
    hold on
    for k = 1:K
        plot(t_start:50, occupancy_array(k,t_start:50), 'LineWidth',2)
    end    
    title(['Occupancy Shares Over Time Stripe ' num2str(bin)])        
    grid on
    legend('Off', 'On1', 'On2')
    saveas(occupancy_fig, [out_dir 'stripe_' num2str(bin) '_' num2str(K) 'state_occupancy_trends.png'], 'png')       
    
    %Calculate 1st order rates from viterbi fits
    rates1 = zeros(K,K);
    for k = 1:K
        dtimes = dwell_times{k};
        s_ct = sum(switch_counts(:,k,:),3);
        ind_vec = 1:K;
        ind_vec = ind_vec(ind_vec~=k);
        s_ct = s_ct(ind_vec);         
        d = -1/mean(dtimes)/dT;
        for n = 1:K-1
            rate = -s_ct(n)/sum(s_ct) * d;
            rates1(ind_vec(n),k) = rate;
        end
        rates1(k,k) = d;   
    end
    stripe_rate_array(:,:,b_ct) = rates1;
    v_struct(b_ct).mean_wt = wt_three;
    v_struct(b_ct).dwell_times = dwell_times;
    v_struct(b_ct).occupancy = occupancy_array;
    v_struct(b_ct).switch_array = switch_counts;
    v_struct(b_ct).R_array = R_array;
    v_struct(b_ct).A_array = A_array;
    v_struct(b_ct).v_array = v_array;
    v_struct(b_ct).binID = bin;
    b_ct = b_ct + 1;
end
save([out_dir num2str(K) 'State_viterbi_dt.mat'],'v_struct');

%%%-----Make Rate Plot-----%%%
colormap((colormap('jet')));
cmap = colormap;
%Make multi-panel fig with fits and originals
rate_fig = figure('Position',[0 0 1024 512]);
increment = floor(60/(2*K));
for j = 1:K
    subplot(1,K,j);
    hold on
    grid on
    index = 1:K;
    index = index(index~=j);            
    legend_cell = {};
    for k = 1:(K-1)
        legend_cell = {legend_cell{:} ['To ' num2str(index(k))]};        
        scatter(stripe_range, reshape(stripe_rate_array(index(k),j,:),1,[]), 50, cmap(increment*((K-1)*(j-1)+k),:), 's', 'filled', 'MarkerEdgeColor', 'black') ;
    end    
    legend(legend_cell{:})
    xlabel('Relative Position')
    ylabel('s^{-1}')
    rate_list = [];
    for k = 1:(K-1)
        r_fit = reshape(stripe_rate_array(index(k),j,:),1,[]);        
        rate_list = [rate_list r_fit];
        p = plot(stripe_range, r_fit,'black','LineWidth',1);        
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    axis([min(stripe_range)-1 max(stripe_range)+1 0 1.25*max(rate_list)]);
    title(['Outflow Rates from State ' num2str(j)]); 
end
saveas(rate_fig, [out_dir 'stripe_' num2str(bin) '_' num2str(K) 'viterbi_rates_by_ap.png'], 'png')       
% saveas(rate_fig, [OutPath '/transition_rates.eps'], 'epsc');
% saveas(rate_fig, [OutPath '/transition_rates.png'], 'png');