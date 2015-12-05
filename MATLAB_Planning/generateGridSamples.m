function [samples] = generateGridSamples(dim, num)

step = 1/num;

space = linspace(step/2,1-step/2,num);
% space = linspace(0,1,num);

xOut = cell(1,dim);
[xOut{:}] = ndgrid(space);

samples = zeros(num^dim,dim);
for i = 1:dim
    xRand = xOut{i} + eps; % need the eps due to some error I haven't diagnosed
    samples(:,i) = xRand(:);
end
end
