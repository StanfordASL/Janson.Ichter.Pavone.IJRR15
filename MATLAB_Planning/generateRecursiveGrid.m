%% Recursive Grid
function [samples] = generateRecursiveGrid(dim,num)
sampleCount = num^dim;

% initialize values
samples = 0.5*ones(1,dim);
count = 1;
step = 1/sqrt(dim^3);
one = ones(1,dim);
% add new axis
while count < sampleCount
    count = count*2;
    step = step/sqrt(dim);
    samples = [bsxfun(@plus,samples,step*one);...
        bsxfun(@minus,samples,step*one)];
end
end