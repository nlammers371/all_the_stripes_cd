%Get eigenvectors & values
eigenvectors = csvread('C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\analysis\shape_modes\eigenvectors.csv');
eigenvalues = csvread('C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\analysis\shape_modes\eigenvalues.csv');
mean_vec = csvread('C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\analysis\shape_modes\mean_vec.csv');


datapath = 'C:\Users\Nicholas\OneDrive\Berkeley Biophysics\projects\all_the_stripes\data\';
load([datapath, 'CompiledParticles_2017_02_08_eve23_4.mat'])

%Get nc14 Traces
fluo_vec_14 = AllTracesVector(:,ncFilter(:,4));
ap_14 = AllTracesAP(ncFilter(:,4));
fluo_vec_14(isnan(fluo_vec_14)) = 0;
fluo_vec_14 = vertcat(fluo_vec_14(219:end,:), zeros(10,size(fluo_vec_14,2)));
%% Conduct Trace Fitting
%Specify number of modes to use (n<=w)
n_modes = 15;
%memory
w = 15;
%Define Traces to be fit
raw_traces = fluo_vec_14(:,5);
%Define fit parameters for each mode
beta_vec = sym('b', [1 n_modes]);
%Filter eigenvector matrix and transpose
shape_mode_mat = eigenvectors(:,end-n_modes + 1:end)';
%Define array to store fit results
shape_fits = struct;
%Define vector to initialize fit paramers
x0 = ones(1,n_modes);
%We will fit each trace piee by piece and then reassemble
for i = 1:(size(raw_traces,2))
    trace_full = raw_traces(:,i)';
    n_fits = length(trace_full) / w;
    out_slice = zeros(n_modes + 2, length(trace_full));
    for j = 1:n_fits
        %Extract appropriate fragment from trace
        vec_full = trace_full((j-1)*w+1:j*w);
        vec_diff = vec_full - vec_full(1);
        %Define objective function
        shape_fit_fun = @(beta_vec)sum((vec_diff - sum(repmat(beta_vec',1,w).* shape_mode_mat) - mean_vec).^2);
        shape_fun = @(beta_vec)mean_vec + sum(repmat(beta_vec',1,w).* shape_mode_mat);
        x = lsqnonlin(shape_fit_fun,x0);
        slice = vertcat(repmat(vec_full(1),1,w),mean_vec,flipud(repmat(x',1,w).*shape_mode_mat));
    	out_slice(:,((j-1)*w+1:j*w)) = slice;
    end
    shape_fits(i).shape_mat = out_slice;
end
%%
for i = 1:15
   
    shape = sum(out_slice(1:12,:));
    plot(1:length(trace_full), trace_full, 1:length(trace_full),shape)
    axis([0 300 0 4.5e5])
    pause(2)

end