function result = ICP(point_clouds_path)
% Input: point_clouds path
% Example: ICP('Data/data/0000000000_normal.pcd')

% Path to the supplemental code
addpath SupplementalCode/

% Read the data
data = readPcd(point_clouds_path);

[ datapoints, dimension ] = size(data);

R = eye(dimension);
t = 0;

result = 0;

end