function visualize(transformed_frames)

figure
for frame_no = 1:length(transformed_frames)
    frame = cell2mat(transformed_frames{frame_no}).';
    scatter3(frame(:, 1),frame(:, 2),frame(:, 3), 0.3), title(string(frame_no))
    hold on
end

end

