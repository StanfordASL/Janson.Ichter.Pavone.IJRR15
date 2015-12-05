%% FMT Function
% Written by Brian Ichter
% Runs the adaptive FMT algorithm and returns the optimal path found
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
    adaptiveFMT(samples, initial, goal, obstacles, r)

global planningType indexMat2Array array2IndexMat gridCount dim;

dim = length(initial);
minR = 0.5*r;
maxR = 1.0*r;
upFactor = 1.05;
downFactor = 0.35;
Cper = 0.7; % adjustment based on collision percentage
Cdist = 1-Cper; % adjustment bsaed on collision distance

% current best, set on completeFMT
% minR = 0.5*r;
% maxR = 1.0*r;
% upFactor = 1.05;
% downFactor = 0.35;

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
    [dist,Nz] = NN(r,z,setV);
end

% remove connections that are not possible
attempts = sum(Nz);
successes = attempts;
for i = 1:length(Nz)
    if Nz(i) == 1
        if(~checkCollision(z,setW(i,:),obstacles))
            Nz(i) = 0;
            successes = successes-1;
        end
    end
end
collisionPercent = 1-successes/attempts;
newR = minR + (r-minR)*(upFactor - collisionPercent*(upFactor-downFactor));

for i = 1:length(Nz)
    if Nz(i) == 1
        insert(tree,nodes(i));
        nodes(i).cost = dist(i);
        nodes(i).r = newR;
    end
end

% keep track of edge costs, add nodes to set H, remove from set W
costEdges = Nz.*dist;
costEdges(isnan(costEdges)) = 0;
costEdges(~costEdges) = NaN;
costAll = costEdges;
setH(Nz,:) = setW(Nz,:);
setW(Nz,:) = NaN;
% go through all nodes
while ~all(z == goal) & any(costEdges)
    [~,loc] = min(costEdges);
    z = setV(loc,:);
    
    % Find points nearby z that have not yet been added
    r = nodes(loc).r;
    if planningType == 3 || planningType == 4 || planningType == 5
        Xnear = indexNNnear(r,loc,z,setW);
    elseif planningType == 2
        [~,Xnear] = NN4(loc,z,setW);
    else
        [~,Xnear] = NN(r,z,setW);
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
            dist = NNdist(r,x,setH);
        end
%         [distSort, idx] = sort(dist+costEdges);
        [~, idx] = min(dist+costEdges);
%         idx(isnan(distSort)) = [];
        
        foundPath = 0;
        newR = minR + (r-minR)*upFactor;
        if newR > maxR
            newR = maxR;
        end
        for j = 1:length(idx)
            costX = costEdges(idx(j)) + dist(idx(j));
            locX = idx(j);
            [valid, freeDist] = checkCollision(x,setH(locX,:),obstacles);
            if(valid)
                if (dist(idx(j)) < r || planningType == 2 || ...
                        planningType == 3 || ...
                        planningType == 4 || planningType == 5)
                    costEdges(i) = costX;
                    costAll(i) = costEdges(i);
                    insert(nodes(locX),nodes(i));
                    nodes(i).cost = costX;
                    foundPath = 1;
                    nodes(i).r = newR;
                end
                break
            else % if collision, downgrade newR
                newR = Cper*(minR + (r-minR)*downFactor)+...
                    Cdist*(newR+freeDist)*0.5;
                nodes(loc).r = newR;
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

    function [distances] = NNdist(~, point, set)
        distances = sqrt(sum((bsxfun(@minus,set,point)).^2,2));
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