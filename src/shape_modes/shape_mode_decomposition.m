% Shape mode bs
datapath = 'C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\data\';
load([datapath, 'CompiledParticles_2017_02_08_eve23_4.mat'])

outpath ='C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\analysis\shape_modes\';

%% Generate Shape n X w Shape Matrix

%Set memory length
w = 15;
%Get nc14 Traces
fluo_vec_14 = AllTracesVector(:,ncFilter(:,4));
ap_14 = AllTracesAP(ncFilter(:,4));
fluo_vec_14(isnan(fluo_vec_14)) = 0;

for i = 1:size(fluo_vec_14,2)
    f_vec = fluo_vec_14(:,i);
    clipped_fluo = f_vec(find(f_vec,1):find(f_vec,1,'last'));
    for j = 1:(length(clipped_fluo) - w + 1)
        shape_raw = clipped_fluo(j:j+w-1);
        shape_diff = shape_raw' - shape_raw(1);
        
        if j == 1 && i == 1
            shape_mat = shape_diff;
            trace_id = i;
        else
            shape_mat = vertcat(shape_mat, shape_diff);
            trace_id = [trace_id i];
        end
    end
end

shape_mat_diffed = shape_mat - repmat(mean(shape_mat),size(shape_mat,1),1);

shape_corr = shape_mat_diffed' * shape_mat_diffed; 

%%
[V, D] = eig(shape_corr);    

csvwrite([outpath 'eigenvectors.csv'], V);
csvwrite([outpath 'eigenvalues.csv'], diag(D));
csvwrite([outpath 'mean_vec.csv'], mean(shape_mat));
% for i = 1:size(V,2)
%     plot(V(:,i))
%     pause(2)
% end
%% Take top 6 shape modes
n_modes = 6;
shape_mode_mat = V(:,end+1-n_modes:end)';
beta_vec = sym('b', [1 n_modes]);
shape_fit_mat = zeros(size(shape_mat));
fit_params = zeros(size(shape_mat,1),n_modes);
avg_shape = mean(shape_mat);
x0 = ones(1,n_modes);
for i = 1:w:w*20
    vec = shape_mat(i,:);
    shape_fit_fun = @(beta_vec)sum((vec - sum(repmat(beta_vec',1,w).* shape_mode_mat) - avg_shape).^2);
    shape_fun = @(beta_vec)avg_shape + sum(repmat(beta_vec',1,w).* shape_mode_mat);
    x = lsqnonlin(shape_fit_fun,x0);
    fit_params(i,:) = x;
    shape_fit_mat(i,:) = shape_fun(x);
end
%%

for i = 1:w:w*20
    plot(1:w,shape_mat(i,:),1:w,shape_fit_mat(i,:));
    pause(1)
end