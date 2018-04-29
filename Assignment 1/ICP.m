function [ opt_trans, full_rms, conv_iter ] = ICP(A1, A2, sampling_method, N_sample, max_iter, show_iter, visualisation)
% ICP               Iterative Closest Point algorithm.
% Input parameters:
% sampling_method   One of the following:
%                   'all': use all data;
%                   'uniform': uniform sampling before main loop;
%                   'random-iter': random sampling at every iteration;
%                   'informative-reg': sampling from informative regions.
% N_sample          Num of points to sample (if sampling_method != 'all').
% max_iter          Maximum iterations to update data (default: 1000).
% show_iter         Boolean for showing training progress (default: false).
% visualisation     Boolean for showing final results (default: true).

% Close all current figures
close all;

% Default parameters
if nargin < 1
    % Read the data
    A1 = load('Data/source.mat');
    A1 = A1.source;
end
if nargin < 2
    A2 = load('Data/target.mat');
    A2 = A2.target;
end
if nargin < 3
   sampling_method = 'informative-reg'; 
end
if nargin < 4
    N_sample = 1000;
end
if nargin < 5
    max_iter = 100;
end
if nargin < 6
   show_iter = false; 
end
if nargin < 7
    visualisation = true;
end

if strcmp(sampling_method, 'all')
    N_sample = size(A1, 2);
end

A1 = A1.';
A1_all = A1; % Used for 'random-iter' sub-sampling
A2 = A2.';
A2_all = A2; % Used for 'random-iter' sub-sampling

% The final transformation
opt_trans = eye(4);

% Check whether N_sample is valid
if N_sample > size(A1, 1)
    N_sample = size(A1, 1);
    disp('N_sample was too large, set to size of A1.')
end
if N_sample > size(A2, 1)
    N_sample = size(A2, 1);
    disp('N_sample was too large, set to size of A2.') 
end

% Uniform sub-sampling
if strcmp(sampling_method, 'uniform')
    ind = randi([1 size(A1, 1)], 1, N_sample);
    A1 = A1(ind, :);
    ind = randi([1 size(A2, 1)], 1, N_sample);
    A2 = A2(ind, :);
end

% Sub-sampling more from informative regions 
if strcmp(sampling_method, 'informative-reg')
    
    [ ~, A1 ] = kmeans(A1, N_sample, 'MaxIter', 500); % Perform k-means
    [ ~, A2 ] = kmeans(A2, N_sample, 'MaxIter', 500); % Perform k-means

    if visualisation
        figure, scatter3(A1(:, 1), A1(:, 2), A1(:, 3), 0.8), title('Informative regions')
        hold on
        scatter3(A2(:, 1), A2(:, 2), A2(:, 3), 0.8)
    end

end



% Calculate RMS of full point cloud
[~, phi] = pdist2(A2_all, A1_all, 'euclidean', 'Smallest', 1);
init_rms = RMS(A1_all, A2_all, phi);

% Visualize both datasets using 3d scatter plot
if visualisation
    figure, scatter3(A1_all(:, 1), A1_all(:, 2), A1_all(:, 3), 0.8), title(strcat({'Init RMS:'},{' '},{num2str(init_rms)}))
    hold on
    scatter3(A2_all(:, 1), A2_all(:, 2), A2_all(:, 3), 0.8)
end

min_rms = 1000;
delta = 0.25;
counter = 0;
% while RMS hasn't converged, update R and t
while counter < max_iter && (counter == 0 || prev_rms <= min_rms + min_rms*delta) % TODO: fix RMS, values are too small ?
    counter = counter + 1;
    if show_iter
        disp("Iteration")
        disp(counter)
    end
    
    % Random sub-sampling in each iteration
    if strcmp(sampling_method, 'random-iter') % TODO: fix, gives poor result
        ind = randi([1 size(A1_all, 1)], 1, N_sample);
        A1 = A1_all(ind, :);
        A1(:, 4) = ones(size(A1, 1), 1);
        A1 = opt_trans * A1.';
        A1 = A1.';
        A1 = A1(:, 1:3);
        ind = randi([1 size(A2_all, 1)], 1, N_sample);
        A2 = A2_all(ind, :);
    end

    % Find the closest point in A2 for each point in A1
    [~, phi] = pdist2(A2, A1, 'euclidean', 'Smallest', 1);

    % Compute the weighted centroids of both point sets:
    w = ones(size(A1, 1), 1); % Set all weights to 1

    p_ = w.' * A1 / sum(w);
    q_ = w.' * A2(phi, :) / sum(w);
   
    % Compute the centered vectors (x = y = 6400x3)
    x = A1 - p_;
    y = A2(phi, :) - q_;
    
    % Compute the d × d covariance matrix (X = Y = 3x6400)
    X = x';
    Y = y';
    W = diag(w);
    S = X * W * Y.';

    % Compute the singular value decomposition
    [U, ~, V] = svd(S);

    % The rotation we are looking for is then
    diag_ones = eye(size(V, 2));
    diag_ones(end, end) = det(V*U.');
    R = V * diag_ones * U.';

    % Compute the optimal translation as
    t = q_.' -  R * p_.';
    
    prev_rms = RMS(A1, A2, phi);
    if prev_rms < min_rms
       min_rms = prev_rms;
       optimal_transformation = opt_trans;
       conv_iter = counter;
    end
    if show_iter
        disp("RMS:")
        disp(prev_rms)
    end
    
    A1 = R * A1.' + t;
    A1 = A1.';
    
    trans = eye(4);
    trans(1:3, 1:3) = R;
    trans(1:3, 4) = t;
    opt_trans = trans * opt_trans;
    
end

A1_all(:, 4) = ones(size(A1_all, 1), 1);
%A1 = opt_trans * A1_all.';
A1 = optimal_transformation * A1_all.';
A1 = A1.';
A2 = A2_all;

% Calculate RMS of full point cloud
A1_all = A1(:, 1:3);
[~, phi] = pdist2(A2_all, A1_all, 'euclidean', 'Smallest', 1);
full_rms = RMS(A1_all, A2_all, phi);
    
if show_iter
    disp('Minimal RMS during covergence:')
    disp(min_rms)
    disp('Amount of iterations until convergence:')
    disp(conv_iter)
    disp('RMS of full point cloud:')
    disp(full_rms)
end

if visualisation
    first = strcat({'Sampling method:'},{' '},{sampling_method}, ...
        {', Num points sampled:'},{' '},{num2str(N_sample)},{'/'},{num2str(size(A1, 1))},{','});
    second = strcat({'RMS:'},{' '},{num2str(full_rms)}, ...
        {', Iterations untill convergence:'},{' '},{num2str(conv_iter)});
    figure, scatter3(A1(:, 1),A1(:, 2),A1(:, 3), 0.8), title({string(first),string(second)});
    hold on
    scatter3(A2(:, 1),A2(:, 2),A2(:, 3), 0.8)
end

end

% Root Mean Square
function rms = RMS(A1, A2, phi)
    rms = sqrt( 1/size(A1, 1) * sum( sqrt(sum(A1 - A2(phi, :)).^2).^2 ) );
end