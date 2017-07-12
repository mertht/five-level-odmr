function Z = radial_gaussian(X, Y, std_dev)
    % Returns a 2-dimensional Gaussian distribution that matches the 
    % dimensions of the X and Y input vectors. The distribution has a mean
    % equal to the center of the X and Y vectors (radius = 0) and a standard
    % deviation of std_dev.
    
    [X2,Y2] = meshgrid(X, Y);
    Z = normpdf(sqrt(X2 .* X2 + Y2 .* Y2), 0, std_dev);
end