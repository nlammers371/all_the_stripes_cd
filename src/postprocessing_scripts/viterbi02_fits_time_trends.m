addpath('../utilities');
%-------------------------- Analyze Viterbi Fits -------------------------%

date_str = '2017-08-16';

%Memory
w = 7;
%State count
ViterbiFolder = ['D:\Data\Nick\projects\all_the_stripes_cd\out\eve7stripes_inf_2017_09_25\2017-09-25\inference_w7'];
out_dir = '..\..\fig\experimental_system\mHMMeve2_orig_inf\2017-08-16\mem8_viterbi_analysis\';
ViterbiFiles = dir(ViterbiFolder);
%If one, single trace plots will be generated
trace_plots = 1;
%If necessary, make output directory
if exist(out_dir) ~= 7
    mkdir(out_dir);
end
trace_dir = [out_dir '/traces/'];
if exist(trace_dir) ~= 7
    mkdir(trace_dir);
end

%list of AP values to check
bin_range = [1:2];
%other inference parameters
graunularity = 5;
dT = 20;
max_time = 60;
b_ct = 1;
%Three State Results
v_fit_struct_three = struct;
for bin = bin_range
    %Track Parameters
    occupancy_array = zeros(3,max_time);
    dwell_times = cell(1,3);
    switch_counts = zeros(3,3,max_time);
    F_counts = cell(1,max_time);
    Fluo_values = cell(1,max_time);
    for i = 1:length(ViterbiFiles)              
        if ~isempty(strfind(ViterbiFiles(i).name,['bin' num2str(bin)]))
            load([ViterbiFolder ViterbiFiles(i).name]);
            %Extract fit parameters
            naive_path = v_fit.z_viterbi;            
            v = v_fit.fluo_viterbi;
            f_exp = v_fit.fluo_exp;
            t_exp = v_fit.time_exp;            
            for m = 1:length(t_exp)
                minute = ceil(t_exp(m)/60 + .01);                                                
                s3 = naive_path(m);                
                f_ct = zeros(1,3);
                f_val = t_exp(m);
                F_t = Fluo_values{minute};
                Fluo_values{minute} = [F_t f_val];
                s_register = naive_path(max(1,m-w+1):m);
                s_register = [s_register zeros(1,w-m-1)];
                for k = 1:K
                    f_ct(k) = sum(s_register == k);
                end
                F_ct = F_counts{minute};
                F_counts{minute} = [F_ct ; f_ct];                
                %Register state occupancy
                occupancy_array(s3,minute) = occupancy_array(s3,minute) + 1;
                %Register switch counts
                if m > 1                    
                    switch_counts(s3_old,s3,minute) = switch_counts(s3_old,s3,minute) + 1;
                end                
                s3_old = s3;
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
                title(['Viterbi Fit for Experimental Trace ' num2str(i) ' (Three States)']);
                grid on
                subplot(2,1,2);
                stairs(naive_path,'Linewidth',1.5)
                grid on
                title(['Inferred Naive Trajectory (' num2str(K) 'States)'])
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
    wt_three = zeros(1,3);
    wait_fig = figure;
    hold on        
    for k = 1:3
        subplot(3,1,k);
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
    occupancy_array = occupancy_array ./ repmat(sum(occupancy_array),3,1);    
    %calculate R(t)
    A_array = zeros(K,K,max_time);
    R_array = zeros(size(A_array));
    window = floor(granularity/2);
    for t = 1:size(occupancy_array,1)
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
    for t = 1:size(occupancy_array,1)
        F = f_counts(t,:);
        fluo = [Fluo_values{max(t-window,1):min(max_time,t+window)}];
        F_ct = vertcat(F_counts{max(t-window,1):min(max_time,t+window)});
        F_square = zeros(K,K);
        b_agg = zeros(1,K);
        for k = 1:K
            F_square(k,:) = sum(F_ct.*repmat(F_ct(:,k),1,K));
            b_agg(k) = sum(fluo'.*F_ct(:,k));
        end
        % solve system                
        v_array(t,:) = linsolve(F_square,b_agg')';        
    end
    occupancy_fig = figure('Position', [0 0 1024 512]);           
    plot(1:max_time, occupancy_array', 'LineWidth',2)
    legend('Off', 'On', 'On2');
    title('Occupancy Shares Over Time')        
    
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
save([out_dir 'ThreeState_viterbi_dt.mat'],'v_struct_three');