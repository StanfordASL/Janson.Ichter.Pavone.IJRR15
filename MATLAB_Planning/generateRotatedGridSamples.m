% generate grid samples and rotate them by a random amount

function [samples] = generateRotatedGridSamples(dim, sideNum)
mult = sqrt(2);
nBound = ceil(sideNum*mult);
multActual = nBound/sideNum;
% numSamples_i = nBound^dim;
allSamples = (generateGridSamples(dim,nBound)-0.5)*multActual;

R = generateRotationMatrix(dim);

% apply rotation
allSamples = (R*(allSamples)')'+0.5;

% apply random offset
% allSamples = bsxfun(@plus,allSamples,rand(1,dim)/sideNum);
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
% numSamples_f = size(allSamples,1);
end