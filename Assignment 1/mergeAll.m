function [] = mergeAll()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

addpath ./SupplementalCode/

% The files
path = './Data/data/';
file = fullfile(path, '*.pcd');
d = dir(file);

file_names = [];
for k = 1:numel(d)
  filename = fullfile(path,d(k).name);
  if ~contains(filename, 'normal')
     file_names = [file_names; filename];
  end
end

for file_no = 1:3%length(file_names)-1
    frame1_filename = file_names(file_no,   :);
    frame2_filename = file_names(file_no+1, :);
    
    frame1 = readPcd(frame1_filename);
    frame1 = frame1(:, 1:3).';
    
    frame2 = readPcd(frame2_filename);
    frame2 = frame2(:, 1:3).';

    merge(frame1, frame2);

end

end

