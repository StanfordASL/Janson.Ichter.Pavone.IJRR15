% generate grid samples and rotate them by a random amount
% takes as input an array of the side lengths (i.e. [4 5 6] for 4 points in
% dim 1, 5 points in dim 2, and 6 points in dim 3

function [samples] = generateWeightedRotatedGridSamples(sides)
dim = size(sides,2);

mult = sqrt(2);
sidesMult = ceil(sides*mult);
multActual = sidesMult./sides;

numSamples_i = prod(sides);
numSamplesMult = prod(sidesMult);

fun = @(A,B) A.*B;
allSamples = bsxfun(fun,generateWeightedGridSamples(sidesMult),multActual);
% offset all points if the sides have enough points
if all(sides > 10)
    allSamples = bsxfun(@plus,allSamples,rand(1,dim)./sides);
end

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

    function ret = givens(n,d1,d2,t)
        ret = eye(n);
        ret([d1 d2],[d1 d2]) = [cos(t) sin(t); -sin(t) cos(t)];
    end
end