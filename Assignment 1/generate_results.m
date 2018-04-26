function [avg_rms, avg_conv_iter] = generate_results(sampling_method, N_runs, ...
    A1, A2, N_sample, max_iter, show_iter, visualisation)

if nargin < 1
    sampling_method = 'uniform';
end
if nargin < 2
    N_runs = 100;
end
if nargin < 3
    % Read the data
    A1 = load('Data/source.mat');
    A1 = A1.source;
end
if nargin < 4
    A2 = load('Data/target.mat');
    A2 = A2.target;
end
if nargin < 5
    N_sample = 6400;
end
if nargin < 6
    max_iter = 100;
end
if nargin < 7
   show_iter = false; 
end
if nargin < 8
    visualisation = false;
end

avg_rms = 0;
avg_conv_iter = 0;
for i = 1:N_runs
    [ ~, min_rms, conv_iter] = ICP(A1, A2, sampling_method, N_sample, max_iter, show_iter, visualisation);
    avg_rms = avg_rms + min_rms;
    avg_conv_iter = avg_conv_iter + conv_iter;
end

avg_rms = avg_rms / N_runs;
avg_conv_iter = avg_conv_iter / N_runs;

disp('The sampling method is:')
disp(sampling_method)
disp('The amount of datapoints used for sampling is:')
disp(N_sample)
disp('The amount of runs using this ICP settings is:')
disp(N_runs)
disp('---RESULTS---')
disp('The average min RMS is:')
disp(avg_rms)
disp('The average convergence iterations is:')
disp(avg_conv_iter)


end
