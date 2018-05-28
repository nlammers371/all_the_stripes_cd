function inf_struct = compile_inference_results(project,inference_type,w,K,Tres,alpha,fluo_type,clipped,...
    clipped_ends,dynamic_bins,t_window,t_inf,DPFolder)

%Compile results into summary structure. Generate summary plots
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = stripe_range(i) - j/3 + 2/3;
    end
end
%Generate filenames and writepath
id_var = [ '/w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/K' num2str(K) '_t_window' num2str(t_window) ...
     '_t_inf' num2str(t_inf) '_' inference_type '/']; 

f_path =  [DPFolder  project  id_var ];
OutPath = ['../../dat/' project  id_var];
FigPath = ['../../fig/experimental_system/' project '/' id_var];
mkdir(OutPath)
mkdir(FigPath)
%---------------------------------Read in Files---------------------------%
files = dir(f_path);
f_names = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,['w' num2str(w)])) && ...
       ~isempty(strfind(files(i).name,['K' num2str(K)]))
        f_names = [f_names {files(i).name}];
    end
end
if isempty(f_names)
    disp(f_path)
    error('No file with specified inference parameters found')
end
%
%Iterate through result sets and concatenate into 1 combined struct
inf_struct = struct;
f_pass = 0;
for f = 1:length(f_names)
    %load the eve validation results into a structure array 'output'        
    load([f_path f_names{f}],'output');        
    if output.skip_flag == 1 
        continue
    elseif output.t_window~=t_window*60 || ~ismember(output.t_inf,t_inf*60)
        continue
    end
    
    for fn = fieldnames(output)'
        inf_struct(f_pass+1).(fn{1}) = output.(fn{1});
    end
    inf_struct(f_pass+1).source = f_names{f};        
    f_pass = f_pass + 1;
    disp(f_pass)
end
% %%
% 
%%% ------------------------------Fig Calculations-------------------------%
%Adjust rates as needed (address negative off-diagonal values)
%Define convenience Arrays and vectors
alpha = inf_struct(1).alpha;
bin_vec = [];
for i = 1:length(inf_struct)
    if length(inf_struct(i).stripe_id) == 21
        bin_vec = [bin_vec 0];
    else        
        bin_vec = [bin_vec mean(inf_struct(i).stripe_id)];
    end
end 
bin_range_vec = unique(bin_vec);
time_vec = [inf_struct.t_inf];
time_index = unique(time_vec);
initiation_rates = zeros(K,length(inf_struct)); % r
pi0_all = zeros(K,length(inf_struct));    
noise_all = zeros(1,length(inf_struct)); % sigma    
dwell_all = NaN(K,length(inf_struct));
% Extract variables and perform rate fitting as needed
for i = 1:length(inf_struct)    
    [initiation_rates(:,i), ranked_r] = sort(60*[inf_struct(i).r]); 
    pi0 = inf_struct(i).pi0(ranked_r);    
    noise_all(i) = sqrt(inf_struct(i).noise); 
    A = reshape(inf_struct(i).A,K,K);
    A = A(ranked_r, ranked_r);
    %Obtain raw R matrix
    R = prob_to_rate(A,Tres);
    inf_struct(i).R_mat = R;
    Rcol = reshape(R,1,[]);
    inf_struct(i).R = Rcol;
    %Check for imaginary and negative elements. If present, perform rate
    %fit   
    inf_struct(i).r_fit_flag = 0;
    if ~isreal(Rcol)
        inf_struct(i).r_fit_flag = 1;
        out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
        inf_struct(i).R_fit = 60*out.R_out;
        r_diag = 60*diag(out.R_out);
    elseif sum(Rcol<0)>K
        inf_struct(i).r_fit_flag = 2;
        R_conv = 60*(R - eye(K).*R);
        R_conv(R_conv<0) = 0;
        R_conv(eye(K)==1) = -sum(R_conv);
        inf_struct(i).R_fit = R_conv;
        r_diag = diag(R_conv);
    else
        inf_struct(i).R_fit = 60*inf_struct(i).R_mat;
        r_diag = 60*diag(inf_struct(i).R_mat); 
    end    
    dwell_all(:,i) = -1./r_diag;
end

% Save rate info
%Make 3D Arrays to stack matrices
R_orig_array = NaN(length(inf_struct),K^2);
R_fit_array = NaN(length(inf_struct),K^2);
A_array = NaN(length(inf_struct),K^2);
for i = 1:length(inf_struct)
    R_fit_array(i,:) = reshape(inf_struct(i).R_fit,1,[]);
    R_orig_array(i,:) = reshape(real(inf_struct(i).R_mat),1,[]);
    A_array(i,:) = reshape(real(inf_struct(i).A_mat),1,[]);
end

% 2D Arrays to store moments (A's are calculated to have for output structure)
avg_R_orig = zeros(length(time_index),K^2,length(bin_range_vec));
std_R_orig = zeros(length(time_index),K^2,length(bin_range_vec));
avg_A = zeros(length(time_index),K^2,length(bin_range_vec));
std_A = zeros(length(time_index),K^2,length(bin_range_vec));
avg_R_fit = zeros(length(time_index),K^2,length(bin_range_vec));
std_R_fit = zeros(length(time_index),K^2,length(bin_range_vec));

for b = 1:length(bin_range_vec)
    bin = bin_range_vec(b);
    for t = 1:length(time_index)
        time = time_index(t);
        A_mean = mean(A_array(bin_vec==bin&time_vec==time,:),1);
        A_reshape = reshape(A_mean,K,K);
        %Normalize A
        A_mean = reshape(A_reshape ./ repmat(sum(A_reshape),K,1),1,[]);
        avg_A(t,:,b) = A_mean;    
        %Calculate and Store R moments
        avg_R_orig(t,:,b) = mean(R_orig_array(bin_vec==bin&time_vec==time,:));
        std_R_orig(t,:,b) = std(R_orig_array(bin_vec==bin&time_vec==time,:));
        avg_R_fit(t,:,b) = mean(R_fit_array(bin_vec==bin&time_vec==time,:));
        std_R_fit(t,:,b) = std(R_fit_array(bin_vec==bin&time_vec==time,:));
    end
end
%Occupancy
occupancy = zeros(K,length(inf_struct));
for i = 1:length(inf_struct)
    [~, ranked_r] = sort([inf_struct(i).r]);
    A = inf_struct(i).A_mat;
    A = A(ranked_r,ranked_r);
    [V,D] = eig(A);
    ind = diag(D)==max(diag(D));
    steady = V(:,ind)./sum(V(:,ind));
    occupancy(:,i) = steady;
    inf_struct(i).occupancy = steady;
end
avg_dwell = NaN(K,length(bin_range_vec),length(time_index));
std_dwell = NaN(K,length(bin_range_vec),length(time_index));
avg_initiation = NaN(K,length(bin_range_vec),length(time_index));
std_initiation = NaN(K,length(bin_range_vec),length(time_index));
avg_occupancy = NaN(K,length(bin_range_vec),length(time_index));
std_occupancy = NaN(K,length(bin_range_vec),length(time_index));
avg_pi0 = NaN(K,length(bin_range_vec),length(time_index));
std_pi0 = NaN(K,length(bin_range_vec),length(time_index));
avg_noise = NaN(1,length(bin_range_vec),length(time_index));
std_noise = NaN(1,length(bin_range_vec),length(time_index));
for i = 1:length(bin_range_vec)
    bin = bin_range_vec(i);
    for t = 1:length(time_index)
        time = time_index(t);
        for k = 1:K
            avg_initiation(k,i,t) = mean(initiation_rates(k,bin_vec==bin&time_vec==time));
            std_initiation(k,i,t) = std(initiation_rates(k,bin_vec==bin&time_vec==time)); 
            avg_occupancy(k,i,t) = mean(occupancy(k,bin_vec==bin&time_vec==time));
            std_occupancy(k,i,t) = std(occupancy(k,bin_vec==bin&time_vec==time)); 
            avg_dwell(k,i,t) = mean(dwell_all(k,bin_vec==bin&time_vec==time));
            std_dwell(k,i,t) = std(dwell_all(k,bin_vec==bin&time_vec==time));        
            avg_pi0(k,i,t) = mean(pi0_all(k,bin_vec==bin&time_vec==time));
            std_pi0(k,i,t) = std(pi0_all(k,bin_vec==bin&time_vec==time));
        end  
        avg_noise(i,t) = mean(noise_all(bin_vec==bin&time_vec==time));
        std_noise(i,t) = std(noise_all(bin_vec==bin&time_vec==time));
    end
end

%%% Make Output Struct With Relevant Fields
hmm_results = struct;
for i = 1:length(bin_range_vec)
    for t = 1:length(time_index)
        ind = (i-1)*length(time_index) + t;
%         hmm_results(ind).N = bin_counts(i);
        hmm_results(ind).initiation_mean = avg_initiation(:,i,t);
        hmm_results(ind).initiation_std = std_initiation(:,i,t);    
        hmm_results(ind).occupancy_mean = avg_occupancy(:,i,t);
        hmm_results(ind).occupancy_std = std_occupancy(:,i,t);        
        hmm_results(ind).pi0_mean = avg_pi0(:,i,t);
        hmm_results(ind).pi0_std = std_pi0(:,i,t);
        hmm_results(ind).dwell_mean = avg_dwell(:,i,t);
        hmm_results(ind).dwell_std = std_dwell(:,i,t);    
        hmm_results(ind).A_mean = avg_A(t,:,i);
        hmm_results(ind).R_orig_mean = avg_R_orig(t,:,i);
        hmm_results(ind).R_orig_std = std_R_orig(t,:,i);
        hmm_results(ind).R_fit_mean = avg_R_fit(t,:,i);
        hmm_results(ind).R_fit_std = std_R_fit(t,:,i);
        hmm_results(ind).noise_mean = avg_noise(i,t);
        hmm_results(ind).noise_std = std_noise(i,t);
        hmm_results(ind).binID = bin_range_vec(i);
        hmm_results(ind).t_inf = time_index(t);
        hmm_results(ind).alpha = alpha;
        hmm_results(ind).dT = Tres;
        hmm_results(ind).fluo_type = fluo_type; 
        hmm_results(ind).clipped = clipped; 
        hmm_results(ind).clipped = clipped_ends; 
    end
end
save([OutPath '/hmm_results_t_window' num2str(t_window)  '_t_inf' num2str(t_inf)  '.mat'],'hmm_results')
disp(['compiled results saved to: ' OutPath])