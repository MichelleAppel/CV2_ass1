function result = ICP(point_clouds_path)
% Input: point_clouds path
% Example: ICP('Data/data/0000000000_normal.pcd')

if nargin < 1
    point_clouds_path = 'Data/data/0000000000_normal.pcd';
end
    
% Path to the supplemental code
addpath SupplementalCode/

% Read the data
data = readPcd(point_clouds_path);

[ datapoints, dimension ] = size(data);

R = eye(dimension);
t = 0;
result = 0;

A1 = load('Data/source.mat');
A1 = A1.source.';
A2 = load('Data/target.mat');
A2 = A2.target.';

% find the closest point in A2 for each point in A1
[~, idx] = pdist2(A1, A2, 'euclidean', 'Smallest', 1);
size(idx);

% visualize both datasets using 3d scatter plot
figure, scatter3(A1(:, 1),A1(:, 2),A1(:, 3))
figure, scatter3(A2(:, 1),A2(:, 2),A2(:, 3))

end