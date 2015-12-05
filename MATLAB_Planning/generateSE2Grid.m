function samples = generateSE2Grid(dim, N, Nt)
mult = sqrt(2);
multN = ceil(N*mult);
multActual = multN/N;

dataRange = (linspace(1/(2*multN), 1 - 1/(2*multN),multN)-0.5)*multActual;
ds = dataRange(2) - dataRange(1);
trange = linspace(0,1,Nt^2);
rotSamples = zeros(N^2*Nt^2,dim);
pt = 1;

R = generateRotationMatrix(dim-1);

for x = dataRange
    for y = dataRange
        for t = 1:length(trange)
            [subx, suby] = ind2sub([Nt,Nt],t);
            subx = subx - (Nt/2+0.5);
            suby = suby - (Nt/2+0.5);
            xy = ([x y] + [subx suby].*[ds ds]/Nt);
            xy = xy*R+0.5;
            rotSamples(pt,:) = [xy(1) xy(2) trange(t)];
            pt = pt+1;
        end
    end
end

goodSamples = true(size(rotSamples,1),1);
for i = length(rotSamples):-1:1
    for j = 1:dim-1
        if(rotSamples(i,j) > 1 || rotSamples(i,j) < 0)
            goodSamples(i) = 0;
            continue;
        end
    end
end
samples = rotSamples(goodSamples,:);
end