function [ transformation ] = merge(frame1, frame2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
addpath ./SupplementalCode

if nargin < 1
    frame1 = readPcd("Data/data/0000000010.pcd");
    frame1 = frame1(:, 1:3).';
end
if nargin < 2
    frame2 = readPcd("Data/data/0000000014.pcd");
    frame2 = frame2(:, 1:3).';
end

% Remove background
frame1 = frame1(:, frame1(3, :) < 2);
frame2 = frame2(:, frame2(3, :) < 2);

[ transformation ] = ICP(frame1, frame2);

end

