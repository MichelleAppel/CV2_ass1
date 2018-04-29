function [ transformation ] = merge(frame1, frame2, sampling_method, N_sample)

addpath ./SupplementalCode

if nargin < 1
    frame1 = readPcd("Data/data/0000000010.pcd");
    frame1 = frame1(:, 1:3).';
end
if nargin < 2
    frame2 = readPcd("Data/data/0000000014.pcd");
    frame2 = frame2(:, 1:3).';
end
if nargin < 3
   sampling_method = 'informative-reg'; 
end
if nargin < 4
    N_sample = 1000;
end

% Remove background
frame1 = frame1(:, frame1(3, :) < 1.42);
frame2 = frame2(:, frame2(3, :) < 1.42);

[ transformation, ~, ~ ] = ICP(frame1, frame2, sampling_method, N_sample);

end

