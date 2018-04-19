function [R, t] = ICP(sampling_method, show_iter, visualisation)
% Input: point_clouds path
% Example: ICP('Data/data/0000000000_normal.pcd')
if nargin < 1
   sampling_method(sampling_method, 'all'); 
end
if nargin < 2
   show_iter = false; 
end
if nargin < 3
    visualisation = false;
end

% Read the data
A1 = load('Data/source.mat');
A1 = A1.source.';
A1 = A1(randperm(size(A1,1)),:);
A2 = load('Data/target.mat');
A2 = A2.target.';
A2 = A2(randperm(size(A2,1)),:);

if strcmp(sampling_method, 'uniform')
    ind = randi([1 size(A1, 1)], 1, 1000);
    A1 = A1(ind, :);
%     ind = randi([1 size(A2, 1)], 1, 1000);
    A2 = A2(ind, :);
end

% Visualize both datasets using 3d scatter plot
if visualisation
    figure, scatter3(A1(:, 1),A1(:, 2),A1(:, 3)), title('A1 START')
    hold on
    scatter3(A2(:, 1),A2(:, 2),A2(:, 3)), title('A2 START')
end

% Init
R = eye(size(A1, 2)); % Initial rotation    of 0
t = 0;                % Initial translation of 0

% Find the closest point in A2 for each point in A1
[~, phi] = pdist2(A2, A1, 'euclidean', 'Smallest', 1);
disp(phi)

prev_rms = 6900000.42; % Intial RMS

counter = 0;
% while RMS hasn't converged, update R and t
while counter < 1
% while RMS(A1, A2, phi) <= prev_rms
    counter = counter + 1;
    if show_iter
        disp("Iteration")
        disp(counter)
    end
    % Find the closest point in A2 for each point in A1
    [~, phi] = pdist2(A2, A1, 'euclidean', 'Smallest', 1);
    
    
    % Compute the weighted centroids of both point sets:
    w = ones(size(A1, 1), 1); % Set all weights to 1
    
    p_ = w.' * A1 / sum(w);
    q_ = w.' * A2 / sum(w);
    
    % Compute the centered vectors (x = y = 6400x3)
    x = A1 - p_;
    y = A2 - q_;
    
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