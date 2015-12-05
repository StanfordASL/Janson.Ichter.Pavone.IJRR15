%% PRM Function
% Written by Brian Ichter
% Runs the PRM algorithm and returns the optimal path found
%
% Given:
%   samples- array of samples (n dimensions by N number of samples)
%   initial- initial state
%   goal- goal state
%   obstacles- obstacle set, given unit length space
%   r- ball radius for NN searches
%
% Return:
%   path- array of the final path from initial to goal state (left empty)
%   tree- array of trees branching from the initial state

function [path, nodes, costAll] = ...
    PRM(samples, initial, goal, obstacles, r, obsType)

% unused output
path = 0;

% create set of all verticies
verts = [samples; initial];
numVerts = size(verts,1);

% create edge lists
edges(1,numVerts) = sledge(-1,-1,NaN);
for i = 1:length(samples)
    edges(i) = sledge(i,-1,NaN);
end

% loop through set of verticies and add edges to list
for i = 1:numVerts
    if(~isnan(verts(i)))
        x = verts(i,:);
        if obsType == 8 || obsType == 9
            [dist, idx] = NNchain(r,x,verts(i+1:numVerts,:));
        else
            [dist, idx] = NN(r,x,verts(i+1:numVerts,:));
        end
        for j = 1:length(idx)
            ji = j+i;
            if(idx(j) ~= 0)
                if(checkCollision(x,verts(ji,:),obstacles,obsType))
                    edgeTo = sledge(i,ji,dist(j));
                    edgeFrom = sledge(ji,i,dist(j));
                    insert(edgeTo,edges(i));
                    insert(edgeFrom,edges(ji));
                end
            end
        end
    end
end

% create nodes
nodes(1,numVerts) = slnode(initial,0,r);
for i = 1:numVerts-1
    nodes(1,i) = slnode(samples(i,:),NaN,r);
end

% Astar search
closedSet(1,numVerts) = nodes(numVerts);
costAll = nan(numVerts,1);
costEst = nan(numVerts,1);

costAll(numVerts) = 0;
costEst(numVerts) = costAll(numVerts) + norm(goal-initial);
currentNode = nodes(numVerts);
while(~all(currentNode.loc == goal) & sum(~isnan(costEst)) > 0)
    [~,currentIdx] = min(costEst); % take lowest estimated cost
    currentNode = nodes(currentIdx);
    currentEdges = edges(currentIdx);
    
    costEst(currentIdx) = NaN;
    closedSet(currentIdx) = currentNode;
    
    while(~isempty(currentEdges.Next)) % loop through all edges
        currentEdges = currentEdges.Next;
        neighborIdx = currentEdges.idxTo;
        neighborNode = nodes(neighborIdx);
        if(isempty(closedSet(neighborIdx).cost)) % if not connecting to closed node
            costTmp = costAll(currentIdx) + currentEdges.cost;
            if(costTmp < costAll(neighborIdx) || isnan(costAll(neighborIdx)))
                insert(currentNode,neighborNode);
                neighborNode.cost = costTmp;
                costAll(neighborIdx) = costTmp;
                costEst(neighborIdx) = costTmp + norm(goal - neighborNode.loc);
            end
        end
    end
end

costAll(numVerts) = []; % return without cost of initial node
% nodes(numVerts) = [];

%% FUNCTION definitions

% nearest neighbors function
    function [distances, nearby] = NN(r, point, set)
        distances = sqrt(sum((bsxfun(@minus,set,point)).^2,2));
        nearby = r > distances;
        distances(~nearby) = nan;
    end

% nearest neighbors function for kinematic chains
    function [distances, nearby] = NNchain(r, point, set)
        diff = abs(bsxfun(@minus,set,point));
        diff_min = bsxfun(@min, diff, 2*pi-diff);
        distances = sum(abs(diff_min),2);
        nearby = r > distances;
        distances(~nearby) = nan;
    end
end