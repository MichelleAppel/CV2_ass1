function RMS = root_mean_square(data1, data2)
    [~, phi] = pdist2(data2, data1, 'euclidean', 'Smallest', 1);   
    RMS = sqrt( 1/size(data1, 1) * sum( sqrt(sum(data1 - data2(phi, :)).^2).^2 ) );
end
