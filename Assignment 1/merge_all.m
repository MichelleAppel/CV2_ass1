function [ transformed_frames ] = merge_all(path, transformations, visualisation)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    path = './Data/data/';
end
if nargin < 2
   t = load('Output/transformations.mat');
   transformations = t.transformations;
end
if nargin < 3
   visualisation = true; 
end

file_names = get_file_names(path);

global_transformation = eye(4);

transformed_frames = cell(length(transformations), 1);

for i = length(transformations):-1:1
    disp(i)
    
    frame1_fileno = cell2mat(transformations{i, 2});

    frame1_filename = file_names(frame1_fileno, :);

    frame = readPcd(frame1_filename);
    frame = frame(:, 1:3).';
    frame = frame(:, frame(3, :) < 2);

    transformation = cell2mat(transformations{i, 1}); % Tranformation
    
    global_transformation = transformation * global_transformation;

    frame(4, :) = ones(size(frame, 2), 1);
    tframe = transformation * frame;
    transformed_frames{i} = num2cell(tframe);
    
end

if visualisation
    visualize(transformed_frames)
end

end

