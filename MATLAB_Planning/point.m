classdef point
    % node A class to represent the contents of a point (node) in the 
    % sample space. The class contains the x and y coordinates, the cost to
    % reach the node, and the ball radius for nearest neighbor searches.
    properties
        x
        y
        cost
        r
    end
    
    methods
        function node = point(x,y,cost,r)
            nargin
            if nargin > 3
                node.x = x;
                node.y = y;
                node.cost = cost;
                node.r = r;
            end
        end
    end
end
        
    