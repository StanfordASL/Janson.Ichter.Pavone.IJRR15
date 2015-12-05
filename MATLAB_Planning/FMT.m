%% FMT Function
% Written by Brian Ichter
% Runs the FMT algorithm and returns the optimal path found
%
% Given:
%   samples- array of samples (n dimensions by N number of samples)
%   initial- initial state
%   goal- goal state
%   obstacles- obstacle set, given unit length space
%   r- ball radius for NN searches
%
% Return:
%   path- array of the final path from initial to goal state
%   tree- array of trees branching from the initial state

function [path, nodes, costAll] = ...
    FMT(samples, initial, goal, obstacles, r, obsType)

global planningType indexMat2Array array2IndexMat gridCount dim;

% set of all verticies excluding the initial point
setV = samples;

% set of nodes not added
setW = samples;

% set of nodes added
setH = NaN(size(samples));

z = initial;
tree = slnode(initial,0,r);

% create nodes
nodes(1,size(samples,1)) = slnode(goal,NaN,r);
for i = 1:size(samples,1)
    nodes(1,i) = slnode(samples(i,:),NaN,r);
end

% set of nearest neighbors
if planningType == 3 || ...
        planningType == 4 || planningType == 5
    [~,initialIdx] = min(sum(((setV - repmat(initial,length(setV),1))).^2,2));
    [dist,Nz] = indexNN(r,initialIdx,initial,setV);
elseif planningType == 2
    [~,initialIdx] = min(sum(((setV - repmat(initial,length(setV),1))).^2,2));
    [dist,Nz] = NN4(initialIdx,initial,setV);
else
    if obsType == 8 || obsType == 9
        [dist,Nz] = NNchain(r,z,setV);
    else
        [dist,Nz] = NN(r,z,setV);
    end
end

% remove connections that are not possible
for i = 1:length(Nz)
    if Nz(i) == 1
        if(~checkCollision(z,setW(i,:),obstacles,obsType)) 
            Nz(i) = 0;
        end
    end
end

for i = 1:length(Nz)
    if Nz(i) == 1
        insert(tree,nodes(i));
        nodes(i).cost = dist(i);
    end
end

% keep track of edge costs, add nodes to set H, remove from set W
% costToGo = sqrt(sum((bsxfun(@minus,setV,goal)).^2,2));
costEdges = Nz.*dist;
costEdges(isnan(costEdges)) = 0;
costEdges(~costEdges) = NaN;
costAll = costEdges;
setH(Nz,:) = setW(Nz,:);
setW(Nz,:) = NaN;
% go through all nodes
tracker = 0; % report number

while ~all(z == goal) && any(costEdges)
    [~,loc] = min(costEdges); % implement as a heap
%     [~,loc] = min(costEdges+costToGo); % with a cost to go heuristic
    z = setV(loc,:);
    
    % Find points nearby z that have not yet been added
    if planningType == 3 || planningType == 4 || planningType == 5
        Xnear = indexNNnear(r,loc,z,setW);
    elseif planningType == 2
        [~,Xnear] = NN4(loc,z,setW);
    else
        if obsType == 8 || obsType == 9
            [~,Xnear] = NNchain(r,z,setW);
        else
            [~,Xnear] = NN(r,z,setW);
        end
    end
    
    indexList = find(Xnear);
    for i_idx = 1:length(indexList)
        i = indexList(i_idx);
        x = setV(i,:);
        
        if planningType == 3 || ...
                planningType == 4 || planningType == 5
            [dist,~] = indexNN(r,i,x,setH);
        elseif planningType == 2
            [dist,~] = NN4(i,x,setH);
        else
            if obsType == 8 || obsType == 9
                [dist,~] = NNchain(r,x,setH);
            else
                [dist,~] = NN(r,x,setH);
            end
        end
        [~, idx] = min(dist+costEdges);
%         [distSort, idx] = sort(dist+costEdges);
%         idx(isnan(distSort)) = [];
        
        foundPath = 0;
        for j = 1:length(idx)
            costX = costEdges(idx(j)) + dist(idx(j));
            locX = idx(j);
            if(checkCollision(x,setH(locX,:),obstacles,obsType))
                if (dist(idx(j)) < r || planningType == 2 || ...
                        planningType == 3 || ...
                        planningType == 4 || planningType == 5)
                    if planningType == 7
                        locPrev = nodes(locX).Prev.loc;
                        if(checkCollision(x,locPrev,obstacles,obsType))
                            % TODO: fix for kinematic chain
                            costX = nodes(locX).Prev.cost + sqrt(sum((locPrev-x).^2));
                            costEdges(i) = costX;
                            costAll(i) = costEdges(i);
                            insert(nodes(locX).Prev,nodes(i));
                            nodes(i).cost = costX;
                            foundPath = 1;
                        else
                            costEdges(i) = costX;
                            costAll(i) = costEdges(i);
                            insert(nodes(locX),nodes(i));
                            nodes(i).cost = costX;
                            foundPath = 1;
                        end
                    else
                        costEdges(i) = costX;
                        costAll(i) = costEdges(i);
                        insert(nodes(locX),nodes(i));
                        nodes(i).cost = costX;
                        foundPath = 1;
                    end
                end
                break
            end
        end
        if ~foundPath
            Xnear(i) = 0;
        end
    end
    
    setH(Xnear,:) = setW(Xnear,:);
    setW(Xnear,:) = NaN;
    
    costEdges(loc) = NaN;
end

path = 0;

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

% nearest neighbors 4 movement lattice function
    function [distances, nearby] = NN4(idx, pos, set)
        %         dim = length(pos);
        %         gridCount = round(length(set)^(1/length(pos)));
        distances = nan(length(set),1);
        nearby = false(length(set),1);
        idxs = array2IndexMat(idx,:)';
        
        offsets = [1 0 -1 0; 0 1 0 -1]';
        for kk = 1:length(offsets)
            ii = idxs(1) + offsets(kk,1);
            jj = idxs(2) + offsets(kk,2);
            if ii > 0 && ii < gridCount + 1 && jj > 0 && jj < gridCount + 1
                newIdx = ii + (jj-1)*gridCount;
                if ~isnan(set(newIdx,1))
                    distances(newIdx) = norm(set(newIdx,:)-pos);
                    nearby(newIdx) = true;
                end
            end
        end
    end

% index nearest neighbors function
% fix to multi dimensions
    function [nearby] = indexNNnear(r, idx, pos, set)
        nearby = false(length(set),1);
        idxDist = round(r*gridCount);
        idxs = array2IndexMat(:,idx);
        
        values = cell(dim,1);
        for ii = 1:dim
            values{ii} = max(idxs(ii)-idxDist,1):min(idxs(ii)+idxDist,gridCount);
        end
        newIdx = indexMat2Array(values{:});
        
        nearby(newIdx(:)) = ~isnan(set(newIdx(:),1));
        
        % keep for down sampling
        distances = nan*set(:,1);
        distances(nearby) = sqrt(sum(bsxfun(@minus,set(nearby,:),pos).^2,2));
        nearby = r > distances;
    end

    function [distances, nearby] = indexNN(r, idx, pos, set)
        %         nearby = indexNNnear(r, idx, pos, set);
        % This part the same as indexNNnear, but faster to split out this way
        
        nearby = false(length(set),1);
        idxDist = round(r*gridCount);
        idxs = array2IndexMat(:,idx);
        
        values = cell(dim,1);
        for ii = 1:dim
            values{ii} = max(idxs(ii)-idxDist,1):min(idxs(ii)+idxDist,gridCount);
        end
        newIdx = indexMat2Array(values{:});
        nearby(newIdx(:)) = ~isnan(set(newIdx(:),1));
        distances = nan*set(:,1);
        distances(nearby) = sqrt(sum(bsxfun(@minus,set(nearby,:),pos).^2,2));
        nearby = r > distances;
        distances(~nearby) = nan;
    end
end
