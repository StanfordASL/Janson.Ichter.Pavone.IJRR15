classdef sledge < handle
    % sledge A class to represent a singly-linked list edge.
    % List edges for a edge via the linked list.
    properties
        idxFrom
        idxTo
        cost
    end
    
    properties(SetAccess = private)
        Next = sledge.empty;
    end
    
    methods
        function edge = sledge(idxFrom,idxTo,cost)
            % Construct a sledge object.
            if nargin > 2
                edge.idxFrom = idxFrom;
                edge.idxTo = idxTo;
                edge.cost = cost;
            end
        end
        
        function insert(newEdge, currentEdge)
            if(~isempty(currentEdge.Next))
                newEdge.Next = currentEdge.Next;
            end
            currentEdge.Next = newEdge;
        end
        
        function removeEdge(edge)
            % Remove a edge from a linked list.
            % unimplemented or used
            if ~isscalar(edge)
                error('Input must be scalar')
            end
        end
        
        function clearList(edge)
            % Clear the list before edge
            next = edge.Next;
            removeEdge(edge)
            while ~isempty(next)
                edge = next;
                next = edge.Next;
                removeEdge(edge)
            end
        end
    end
    
    methods (Access = private)
        function delete(edge)
            % Delete all edges
%             clearList(edge)
        end
    end % private methods
end % classdef