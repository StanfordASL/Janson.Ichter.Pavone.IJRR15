classdef slnode < handle
    % slnode A class to represent a singly-linked list node.
    % Link multiple slnode objects together to create linked lists.
    properties
        loc
        cost
        r
    end
    properties(SetAccess = private)
        Prev = slnode.empty;
    end
    
    methods
        function node = slnode(loc,cost,r)
            % Construct a slnode object.
            if nargin > 2
                node.loc = loc;
                node.cost = cost;
                node.r = r;
            end
        end
        
        function insert(newNode, nodeAfter)
            nodeAfter.Prev = newNode;
        end
        
        function removeNode(node)
            % Remove a node from a linked list.
            % unimplemented or used
            if ~isscalar(node)
                error('Input must be scalar')
            end
        end
        
        function [solutionPath] = smoothPath(solutionPath, obstacles, obsType)
            for i = 1:2
                node1 = solutionPath;
                nodeIntermediate = node1.Prev;
                while(~isempty(nodeIntermediate) && ~isempty(nodeIntermediate.Prev))
                    node2 = nodeIntermediate.Prev;
                    if checkCollision(node1.loc, node2.loc, obstacles, obsType);
                        insert(node2,node1);
                    end
                    node1 = node1.Prev;
                    nodeIntermediate = node1.Prev;
                end
                node1 = solutionPath;
                if(~isempty(node1) && ~isempty(node1.Prev))
                    node2 = node1.Prev;
                    node3 = node2.Prev;
                    while(~isempty(node1) && ~isempty(node1.Prev) ...
                            && ~isempty(node1.Prev.Prev))
                        
                        shortCutStart = rand();
                        shortCutEnd = rand();
                        
                        midnode12loc = (1-shortCutStart)*node1.loc + ...
                            shortCutStart*node2.loc;
                        midnode23loc = (1-shortCutEnd)*node2.loc + ...
                            (shortCutEnd)*node3.loc;
                        
                        if checkCollision(midnode12loc, midnode23loc, obstacles, obsType);
                            midnode12 = slnode(midnode12loc,1,1);
                            midnode23 = slnode(midnode23loc,1,1);
                            
                            node1.Prev = midnode12;
                            midnode12.Prev = midnode23;
                            midnode23.Prev = node3;
                            
                            node1 = node3;
                        else
                            node1 = node2;
                        end
                        
                        if(~isempty(node1.Prev) && ~isempty(node1.Prev.Prev))
                            node2 = node1.Prev;
                            node3 = node2.Prev;
                        end
                    end
                end
            end
        end
        
        function updateCost(node)
            node.cost = 0;
            currentNode = node;
            while ~isempty(currentNode.Prev)
                node.cost = node.cost + ...
                    norm(currentNode.loc - currentNode.Prev.loc);
                currentNode = currentNode.Prev;
            end
        end
        
        function clearList(node)
            % Clear the list before node
            prev = node.Prev;
            removeNode(node)
            while ~isempty(prev)
                node = prev;
                prev = node.Prev;
                removeNode(node)
            end
        end
    end
    
    methods (Access = private)
        function delete(node)
            % Delete all nodes
            clearList(node)
        end
    end % private methods
end % classdef