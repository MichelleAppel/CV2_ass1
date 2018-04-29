function [ transformations ] = comp_transformations(path, step_size, ...
    print_step, sampling_method, N_sample)

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
   sampling_method = 'uniform'; 
end
if nargin < 5
    N_sample = 1000;
end

file_names = get_file_names(path);

transformations = cell(length(1:step_size:length(file_names)-step_size), 3);

for file_no = 1:step_size:length(file_names)-step_size
    if print_step
        fprintf(strcat("\nStep:", string((file_no-1)/step_size + 1), "/", string(length(1:step_size:length(file_names)-1))))
    end
    
    frame1_filename = file_names(file_no,           :);
    frame2_filename = file_names(file_no+step_size, :);
    
    frame1 = readPcd(frame1_filename);
    frame1 = frame1(:, 1:3).';
    
    frame2 = readPcd(frame2_filename);
    frame2 = frame2(:, 1:3).';

    [ transformation ] = merge(frame1, frame2, sampling_method, N_sample);
    transformations{(file_no-1)/step_size + 1, 1} = num2cell(transformation);
    transformations{(file_no-1)/step_size + 1, 2} = num2cell(file_no);
    transformations{(file_no-1)/step_size + 1, 3} = num2cell(file_no+step_size);
end

file_name = strcat('Output/transformations_step_', num2str(step_size), ...
    '_', sampling_method, '_N_', num2str(N_sample), '.mat');
save(file_name, 'transformations');

end

