function [opt_R, opt_t] = ICP(sampling_method, N_sample, max_iter, show_iter, visualisation)
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
%%% Input: point_clouds path
%%% Example: ICP('Data/data/0000000000_normal.pcd')

close all;
clear;
clc;

% Default parameters
if nargin < 1
   sampling_method = 'all'; 
end
if nargin < 2
    N_sample = 2000;
end
if nargin < 3
    max_iter = 10;
end
if nargin < 4
   show_iter = false; 
end
if nargin < 5
    visualisation = true;
end

% Read the data
A1 = load('Data/source.mat');
A1 = A1.source.';
%A1_all = A1; % Used for 'random-iter' sub-sampling
A2 = load('Data/target.mat');
A2 = A2.target.';
%A2_all = A2; % Used for 'random-iter' sub-sampling

opt_R = eye(3);
opt_t = [0; 0; 0];

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
    % TODO
end

% Visualize both datasets using 3d scatter plot
if visualisation
    figure, scatter3(A1(:, 1), A1(:, 2), A1(:, 3)), title('A1 START')
    hold on
    scatter3(A2(:, 1), A2(:, 2), A2(:, 3)), title('A2 START')
end

counter = 0;
% while RMS hasn't converged, update R and t
while counter < max_iter    % while RMS(A1, A2, phi) <= prev_rms
    counter = counter + 1;
    if show_iter
        disp("Iteration")
        disp(counter)
    end
    
    % Random sub-sampling in each iteration
    %if strcmp(sampling_method, 'random-iter')
    %    ind = randi([1 size(A1_all, 1)], 1, N_sample);
    %    A1 = A1_all(ind, :);
    %    ind = randi([1 size(A2_all, 1)], 1, N_sample);
    %    A2 = A2_all(ind, :);
    %end

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
    if show_iter
        disp("RMS:")
        disp(prev_rms)
    end
    
    A1 = R * A1.' + t;
    A1 = A1.';
    
    opt_R = opt_R * R;
    opt_t = opt_t + t;
  
end

if visualisation
    figure, scatter3(A1(:, 1),A1(:, 2),A1(:, 3)), title('A1 END')
    hold on
    scatter3(A2(:, 1),A2(:, 2),A2(:, 3)), title('A2 END')
end
    
end

% Root Mean Square
function rms = RMS(A1, A2, phi)
    rms = sqrt( 1/size(A1, 1) * sum( sqrt(sum(A1 - A2(phi, :)).^2).^2 ) );
end