function [samples] = generateHexGridSamples(dim, sampleCount)

if(dim ~= 2)
    disp('To use generateHexGridSamples, must be 2 dimensional');
    return
end

numX = round(sqrt(sampleCount*sqrt(3)/2));
numY = round(sqrt(sampleCount*2/sqrt(3)));

stepX = 1/numX;
stepY = 1/numY;
spaceX = linspace(stepX/2,1-stepX/2,numX);
spaceY = linspace(stepY/2,1-stepY/2,numY);

xOut = cell(1,dim);
[xOut{:}] = meshgrid(spaceX,spaceY);

for i = 1:2:numY
    xOut{1}(i,:) = xOut{1}(i,:) + stepX/4;
end
for i = 2:2:numY
    xOut{1}(i,:) = xOut{1}(i,:) - stepX/4;
end

samples = zeros(numX*numY,dim);
for i = 1:dim
    xRand = xOut{i} + eps;
    samples(:,i) = xRand(:);
end
end
