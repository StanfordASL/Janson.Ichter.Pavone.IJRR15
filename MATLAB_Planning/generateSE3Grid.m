% generate grid samples and rotate them by a random amount
% takes as input an array of the side lengths (i.e. [4 5 6] for 4 points in
% dim 1, 5 points in dim 2, and 6 points in dim 3

function [samples] = generateSE3Grid(sides)
dim = length(sides);

mult = sqrt(2);
sidesMult = ceil(sides*mult);
multActual = sidesMult./sides;

numSamples_i = prod(sides);
numSamplesMult = prod(sidesMult);

fun = @(A,B) A.*B;
allSamples = bsxfun(fun,generateWeightedGridSamples(sidesMult),multActual);
allSamples = bsxfun(@plus,allSamples,rand(1,dim)./sidesMult); % offset
allSamples(:,1) = allSamples(:,1) + 1/(sidesMult(1))*(allSamples(:,4)-0.5); % split out dimensions
allSamples(:,2) = allSamples(:,2) + 1/(sidesMult(2))*(allSamples(:,5)-0.5);
allSamples(:,3) = allSamples(:,3) + 1/(sidesMult(3))*(allSamples(:,6)-0.5);

R = generateRotationMatrix(dim);

center = multActual/2;
allSamples = bsxfun(@plus,(R*bsxfun(@minus,allSamples,center)')',0.5);

validSamples = true(size(allSamples,1),1);
for i = size(allSamples,1):-1:1
    for d = 1:dim
        if(allSamples(i,d) > 1 || allSamples(i,d) < 0)
            validSamples(i) = false;
            continue;
        end
    end
end
samples = allSamples(validSamples,:);
numSamples_f = size(allSamples,1);
end