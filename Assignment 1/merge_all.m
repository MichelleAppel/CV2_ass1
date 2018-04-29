function [ transformed_frames ] = merge_all(path, transformations, visualisation)

addpath ./SupplementalCode/

if nargin < 1
    path = './Data/data/';
end
if nargin < 2
   t = load('Output/transformations_step_10_uniform_N_100');
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
    tframe = global_transformation * frame;
    transformed_frames{i} = num2cell(tframe);
end

avg_rms = 0;
for i=2:length(transformed_frames)
    frame1 = cell2mat(transformed_frames{i-1});
    frame2 = cell2mat(transformed_frames{i});
    frame1 = frame1(1:3, :).';
    frame2 = frame2(1:3, :).';
    rms = root_mean_square(frame1, frame2);
    avg_rms = avg_rms + rms;
    disp(rms)
end
avg_rms = avg_rms / length(transformations);
disp('Average RMS:')
disp(avg_rms)

if visualisation
    visualize(transformed_frames)
end

end

