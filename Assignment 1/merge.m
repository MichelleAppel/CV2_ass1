function [] = merge(frame1, frame2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
addpath ./SupplementalCode

%pointCloud = readPcd('Data/data/0000000099.pcd');
%disp(pointCloud)

if nargin < 1
    frame1 = readPcd("Data/data/0000000012.pcd");
    frame1 = frame1(:, 1:3).';
    size(frame1)
end
if nargin < 2
    frame2 = readPcd("Data/data/0000000014.pcd");
    frame2 = frame2(:, 1:3).';
    size(frame2)
end

% Remove background
frame1 = frame1(:, frame1(3, :) < 2);
frame2 = frame2(:, frame2(3, :) < 2);

% Plot
%frame1 = frame1.';
%frame2 = frame2.';
%figure, scatter3(frame1(:, 1),frame1(:, 2),frame1(:, 3)), title('frame1 END')
%hold on
%scatter3(frame2(:, 1),frame2(:, 2),frame2(:, 3)), title('frame2 END')

ICP(frame1, frame2)


end

