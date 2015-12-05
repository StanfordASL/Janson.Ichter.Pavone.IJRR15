function [pos] = chainPosition(thetas)
    dim = size(thetas,2);
    linkLength = 1/dim/sqrt(2);
    theta = 0;
    w = ones(1,2)*0.5;
    for link = 1:dim
        theta = theta + thetas(link);
        v = w;
        w = v + [cos(theta) sin(theta)]*linkLength;
    end
    pos = w;
end
