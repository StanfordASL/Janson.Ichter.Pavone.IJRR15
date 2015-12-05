function [samples] = generatePertGridSamples(dim, num)

step = 1/num;

space = linspace(step/2,1-step/2,num);
% space = linspace(0,1,num);

xOut = cell(1,dim);
[xOut{:}] = ndgrid(space);

samples = zeros(num^dim,dim);
for i = 1:dim
    xRand = xOut{i} + rand(size(xOut{i}))*step-step/2;
    samples(:,i) = xRand(:);
end
end
