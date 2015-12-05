classdef treeNode < handle
    % treeNode A class to represent a tree node. 
    
    properties
        loc
        cost
        r
        index
        parent = treeNode.empty;
        children = treeNode.empty;
    end
    
    methods
        function node = treeNode(loc,cost,r,index)
            % constructor
            if nargin > 3
                node.loc = loc;
                node.cost = cost;
                node.r = r;
                node.index = index;        
            end
        end
        
        function insertNode(parentNode, childNode)
            parentNode.children(length(parentNode.children)+1) = childNode;
            childNode.parent = parentNode;
        end
    end
end