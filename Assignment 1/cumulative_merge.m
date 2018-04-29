function [ complete_point_cloud ] = cumulative_merge(path, step_size, print_step, visualize)

addpath ./SupplementalCode/

if nargin < 1
    path = './Data/data/';
end
if nargin < 2
    step_size = 1;
end
if nargin < 3
   print_step = true; 
end
if nargin < 4
   visualize = false; 
end

file_names = get_file_names(path);

frame1 = readPcd(file_names(1, :));
frame1 = frame1(:, 1:3).';
complete_point_cloud = frame1;

for file_no = 1:step_size:length(file_names)-step_size
    if print_step
        fprintf(strcat("\nStep:", string((file_no-1)/step_size + 1), "/", string(length(1:step_size:length(file_names)-1))))
    end
    
    frame1 = complete_point_cloud;
    
    frame2_filename = file_names(file_no+step_size, :);
    frame2 = readPcd(frame2_filename);
    frame2 = frame2(:, 1:3).';

    [ transformation ] = merge(frame1, frame2);
   
    frame1(4, :) = ones(size(frame1, 2), 1);
    tframe1 = transformation * frame1;
    tframe1 = tframe1(1:3, :);
    
    complete_point_cloud = [tframe1, frame2];
    complete_point_cloud = complete_point_cloud(:, complete_point_cloud(3, :) < 1.42);
    
    
%     complete_point_cloud(:, end:end+size(tframe2, 2)) = tframe2;
%     complete_point_cloud = [tframe1, frame2];
    
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

if visualize
    complete_point_cloud = complete_point_cloud.';
    scatter3(complete_point_cloud(:, 1),complete_point_cloud(:, 2),complete_point_cloud(:, 3), 0.3)
end

end

