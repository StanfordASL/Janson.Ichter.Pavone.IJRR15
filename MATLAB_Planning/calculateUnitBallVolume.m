function [vol] = calculateUnitBallVolume(dim)
if (dim==0)
    vol = 1;
    return
elseif (dim==1)
    vol = 2;
    return
end
vol = 2*pi/dim*calculateUnitBallVolume(dim-2);
end