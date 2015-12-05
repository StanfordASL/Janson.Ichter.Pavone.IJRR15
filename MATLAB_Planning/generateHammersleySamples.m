function [samples] = generateHammersleySamples(dim, sampleCount)
% Generate samples using a Hammersley Sequence

primeNums = primes(1000);
if dim > length(primeNums)
    disp('Change the max prime number, in generateHaltonSamples()');
end
base = primeNums(1:dim-1);

samples = zeros(sampleCount,dim);

samples(:,1) = (1:sampleCount)'/sampleCount;
for idx=1:sampleCount
    k = idx + 20;
    p = base;
    while k > 0
        a = mod(k,base);
        samples(idx,2:dim) = samples(idx,2:dim) + a./p;
        k = floor(k./base);
        p = p.*base;
    end
end
samples = samples+eps;

end

% plot3(samples(:,1),samples(:,2),samples(:,3),'.')
% plot(samples(:,1),samples(:,2),'.')