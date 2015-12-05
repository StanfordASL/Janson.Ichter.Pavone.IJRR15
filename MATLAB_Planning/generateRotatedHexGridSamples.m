function [samples] = generateRotatedHexGridSamples(dim, sampleCount)

if(dim ~= 2)
    disp('To use generateHexGridSamples, must be 2 dimensional');
    return
end

mult = 2;
nBound = ceil(sampleCount*mult);
multActual = nBound/sampleCount;
% numSamples_i = nBound^dim
allSamples = (generateHexGridSamples(dim,nBound)-0.5)*sqrt(multActual);
allSamples = bsxfun(@plus,allSamples,rand(1,dim)/sqrt(sampleCount));

R = generateRotationMatrix(dim);

allSamples = (R*(allSamples)')'+0.5;
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
end
