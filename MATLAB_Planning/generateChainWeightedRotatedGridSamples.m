% generate grid samples and rotate them by a random amount
% takes as input an array of the side lengths (i.e. [4 5 6] for 4 points in
% dim 1, 5 points in dim 2, and 6 points in dim 3

function [samples] = generateChainWeightedRotatedGridSamples(sides)
dim = size(sides,2);
samples = generateWeightedGridSamples(sides);
R = generateRotationMatrix(dim);
samples = mod((R*samples')',1);
end