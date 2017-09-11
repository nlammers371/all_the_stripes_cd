addpath('../utilities/');

date_str = '2017_06_15';
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
all_outputs = struct;
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'
    load(['../../inference_results/' date_str '/' filenames{f}], 'output');
    all_outputs(f).results = output;
end


% ---------- Sort the states in an ascending emission order ------------
for i = 1:length(all_outputs)
    r_temp = all_outputs(i).results.r;
    [r_sort, I] = sort(r_temp);

    all_outputs(i).results.r = r_sort;

    v_temp = all_outputs(i).results.v;
    all_outputs(i).results.v = v_temp(I);

    R_temp = all_outputs(i).results.R_mat;
    R_temp = R_temp(I,I);
    all_outputs(i).results.R_mat = R_temp;
    all_outputs(i).results.R = R_temp(:);

    A_temp = all_outputs(i).results.A_mat;
    A_temp = A_temp(I,I);
    all_outputs(i).results.A_mat = A_temp;
    all_outputs(i).results.A = A_temp(:);

    A_log_temp = all_outputs(i).results.A_mat;
    A_log_temp = A_log_temp(I,I);
    all_outputs(i).results.A_log = A_log_temp;

    pi0_temp = all_outputs(i).results.pi0;
    all_outputs(i).results.pi0 = pi0_temp(I);

    pi0_log_temp = all_outputs(i).results.pi0_log;
    all_outputs(i).results.pi0_log = pi0_log_temp(I);
end
% -------------------- Best Inferred Params set parameters --------------


% number of data sets per trace count
outputs_best = [];
for i = 1:length(all_outputs)
    out_best = [all_outputs(i).AP, all_outputs(i).results.A', all_outputs(i).results.R', all_outputs(i).results.v', all_outputs(i).results.noise, all_outputs(i).results.pi0,   ...
             all_outputs(i).results.N all_outputs(i).results.total_steps, all_outputs(i).results.total_time];

    outputs_best = vertcat(real(out_best), outputs_best);  
end

csvwrite([outpath  'best_results.csv'], outputs_best);