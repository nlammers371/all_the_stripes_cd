addpath('../utilities');
%-------------------------- Analyze Viterbi Fits -------------------------%

date_str = '2017-08-16';

%Memory
w = 8;
%State count
ThreeStateFolder = '..\..\fig\experimental_system\mHMMeve2_orig_inf\2017-08-16\mem8_states3_s0_35\viterbi_fits\';
TwoStateFolder = '..\..\fig\experimental_system\mHMMeve2_orig_inf\2017-08-16\mem8_states2_s0_35\viterbi_fits\';
out_dir = '..\..\fig\experimental_system\mHMMeve2_orig_inf\2017-08-16\mem8_viterbi_analysis\';
TwoFiles = dir(TwoStateFolder);
ThreeFiles = dir(ThreeStateFolder);
if length(TwoFiles) ~= length(ThreeFiles)
    warning('Inconsistent Trace Counts Between Comparison Architectures')
end
%If one, single trace plots will be generated
trace_plots = 0;
%If necessary, make output directory
if exist(out_dir) ~= 7
    mkdir(out_dir);
end
trace_dir = [out_dir '/traces/'];
if exist(trace_dir) ~= 7
    mkdir(trace_dir);
end
%list of AP values to check
bin_list = [-1:1];
dT = 22.5;
max_time = 35;
b_ct = 1;
%Struct to store two state results
v_fit_struct_two = struct;
%Three State Results
v_fit_struct_three = struct;
for bin = -1:1
    %Track wait times
    occupancy_vec_three = zeros(3,max_time);
    dwell_times_three = cell(1,3);
    switch_counts_three = zeros(3,3,max_time);
    
    occupancy_vec_two = zeros(2,max_time);
    dwell_times_two = cell(1,2);
    switch_counts_two = zeros(2,2,max_time);
    
    for i = 1:length(TwoFiles)
        if isempty(strfind(TwoFiles(i).name,['bin' num2str(bin)])) ~= 1
            %load in data structure named "v_fit" 
            load([TwoStateFolder TwoFiles(i).name]);
            %Make plot comparison fig
            two_naive_path = v_fit.z_viterbi;            
            two_v = v_fit.fluo_viterbi;
            two_exp = v_fit.fluo_exp;
            load([ThreeStateFolder ThreeFiles(i).name]);
            %Make plot comparison fig
            three_naive_path = v_fit.z_viterbi;            
            three_v = v_fit.fluo_viterbi;
            three_exp = v_fit.fluo_exp;
            time_exp = v_fit.time_exp;
            for m = 1:length(time_exp)
                minute = ceil(time_exp(m)/60);
                s2 = two_naive_path(m);
                occupancy_vec_two(s2,minute) = occupancy_vec_two(s2,minute) + 1;
                
                s3 = three_naive_path(m);
                occupancy_vec_three(s3,minute) = occupancy_vec_three(s3,minute) + 1;
                if m > 1
                    switch_counts_two(s2_old,s2,minute) = switch_counts_two(s2_old,s2,minute) + 1;
                    switch_counts_three(s3_old,s3,minute) = switch_counts_three(s3_old,s3,minute) + 1;
                end
                s2_old = s2;
                s3_old = s3;
            end
            if trace_plots
                cp_fig = figure('Position', [0 0 1024 1024], 'Visible', 'off');            
                subplot(2,2,1);
                hold on
                plot(three_exp, 'Linewidth', 1.5);
                plot(three_v, 'Linewidth', 1.5);
                legend('Raw Trace', 'Viterbi Fit');
                ylabel('AU')
                xlabel('Time Steps');
                title(['Viterbi Fit for Experimental Trace ' num2str(i) ' (Three States)']);
                grid on
                subplot(2,2,2);
                hold on
                plot(two_exp, 'Linewidth', 1.5);
                plot(two_v, 'Linewidth', 1.5);
                legend('Raw Trace', 'Viterbi Fit');
                ylabel('AU')
                xlabel('Time Steps');
                title(['Viterbi Fit for Experimental Trace  ' num2str(i) ' (Two States)']);
                grid on
                subplot(2,2,3);
                stairs(three_naive_path,'Linewidth',1.5)
                grid on
                title('Inferred Naive Trajectory (Three States)')

                subplot(2,2,4);
                stairs(two_naive_path,'Linewidth',1.5)
                grid on
                title('Inferred Naive Trajectory (Two States)')

                saveas(cp_fig, [trace_dir 'bin_' num2str(bin) '_trace_' num2str(i) '.png'], 'png')            
            end
            %Add dwell times for each state. Ignore trace boundaries
            for k = 1:length(dwell_times_two)
                d_times = dwell_times_two{k};
                binary = two_naive_path == k;
                d_binary = diff(binary);
                starts = find(d_binary == 1) + 1;
                stops = find(d_binary == -1);
                stops = stops(stops>=min(starts));
                starts = starts(starts<=max(stops));
                gaps = stops - starts + 1;
                dwell_times_two{k} = [d_times gaps];
            end               
            for k = 1:length(dwell_times_three)
                d_times = dwell_times_three{k};
                binary = three_naive_path == k;
                d_binary = diff(binary);
                starts = find(d_binary == 1) + 1;
                stops = find(d_binary == -1);
                stops = stops(stops>=min(starts));
                starts = starts(starts<=max(stops));
                gaps = stops - starts + 1;
                dwell_times_three{k} = [d_times gaps];
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
        wt_ct = histc(dwell_times_three{k}, wait_bins);        
        wt_mean = sum(wait_bins.*wt_ct)/sum(wt_ct);
        wt_three(k) = wt_mean;
        wt_ct = wt_ct / sum(wt_ct);
        fit = poisspdf(wait_bins,wt_mean);
        bar(wait_bins*dT, wt_ct,'BarWidth', 1,'FaceAlpha',0.5);
        plot(wait_bins*dT,fit,'LineWidth',2)        
        title(strvcat(['Distribution of Dwell Times with Poisson Fit: State ' num2str(k-1)],...
                    ['']))             
        grid on        
    end
    saveas(wait_fig, [out_dir 'bin_' num2str(bin) '_3state_fit.png'], 'png')   
    
    wt_two = zeros(1,2);
    wait_fig = figure;
    hold on        
    for k = 1:2
        subplot(2,1,k);
        hold on
        wait_bins = 1:ceil(5*60/dT);
        wt_ct = histc(dwell_times_two{k}, wait_bins);  
        wt_mean = sum(wait_bins.*wt_ct)/sum(wt_ct);
        wt_two(k) = wt_mean;
        wt_ct = wt_ct / sum(wt_ct);
        fit = poisspdf(wait_bins,wt_mean);
        bar(wait_bins*dT, wt_ct,'BarWidth', 1,'FaceAlpha',0.5);
        plot(wait_bins*dT,fit,'LineWidth',2)
%         f = histfit(wt_ct,length(wait_bins),'poisson');
        title(strvcat(['Distribution of Dwell Times with Poisson Fit: State ' num2str(k-1)],...
                    ['']))             
        grid on        
    end
    saveas(wait_fig, [out_dir 'bin_' num2str(bin) '_2state_fit.png'], 'png')
    
    occupancy_fig = figure('Position', [0 0 1024 512]);
    subplot(2,1,1);    
    o3_normed = occupancy_vec_three ./ repmat(sum(occupancy_vec_three),3,1);
    plot(1:max_time, o3_normed', 'LineWidth',2)
    legend('Off', 'On', 'On2');
    title('Occupancy Shares Over Time')
    subplot(2,1,2);
    o2_normed = occupancy_vec_two ./ repmat(sum(occupancy_vec_two),2,1);
    plot(1:max_time, o2_normed', 'LineWidth',2)
    legend('Off', 'On');
    saveas(occupancy_fig, [out_dir 'bin_' num2str(bin) '_occupancy_trends.png'], 'png')
    v_struct_two(b_ct).mean_wt = wt_two;
    v_struct_two(b_ct).dwell_times = dwell_times_two;
    v_struct_two(b_ct).occupancy = occupancy_vec_two;
    v_struct_two(b_ct).switch_array = switch_counts_two;
    
    v_struct_three(b_ct).mean_wt = wt_three;
    v_struct_three(b_ct).dwell_times = dwell_times_three;
    v_struct_three(b_ct).occupancy = occupancy_vec_three;
    v_struct_three(b_ct).switch_array = switch_counts_three;
    b_ct = b_ct + 1;
end
save([out_dir 'bin_' num2str(bin) '_2state_viterbi_dt.mat'],'v_struct_two');
save([out_dir 'bin_' num2str(bin) '_3state_viterbi_dt.mat'],'v_struct_three');