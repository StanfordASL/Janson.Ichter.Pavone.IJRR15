% generate grid samples and rotate them by an arbitrary amount
% takes as input an array of the side lengths (i.e. [4 5 6] for 4 points in
% dim 1, 5 points in dim 2, and 6 points in dim 3

function [samples] = generateWeightedGridSamples(sides)

dim = length(sides);

spaces = cell(1,dim);
num = 1;
for i = 1:dim
    sideNum = sides(i);
    step = 1/sideNum;
    spaces{i} = linspace(step/2,1-step/2,sideNum);
    num = num*sideNum;
end

xOut = cell(1,dim);
[xOut{:}] = ndgrid(spaces{:});
samples = zeros(num,dim);
for i = 1:dim
    xRand = xOut{i} + eps; % need eps due to a numerical error
    samples(:,i) = xRand(:);
end
end
