function [ file_names ] = get_file_names(path)

file = fullfile(path, '*.pcd');
d = dir(file);

file_names = [];
for k = 1:numel(d)
  filename = fullfile(path,d(k).name);
  if ~contains(filename, 'normal')
     file_names = [file_names; filename];
  end
end


end

