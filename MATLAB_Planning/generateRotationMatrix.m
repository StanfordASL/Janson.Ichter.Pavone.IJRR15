function [R] = generateRotationMatrix(dim);

% random rotation
% [R,~] = qr(rand(dim));

% rotation to [1 1 .. 1]
% y = ones(1,dim)';
% v = y/norm(y);
% [q,r] = qr(v);
% R = q*r(1);

theta = 10*pi*pi/180;
R = eye(dim);
if dim > 2
    for i = 1:dim-1
        R = R*givens(dim,i,i+1,theta);
    end
end
R = R*givens(dim,dim,1,theta);
% R = R*R;

    function ret = givens(n,d1,d2,t)
        ret = eye(n);
        ret([d1 d2],[d1 d2]) = [cos(t) sin(t); -sin(t) cos(t)];
    end

end