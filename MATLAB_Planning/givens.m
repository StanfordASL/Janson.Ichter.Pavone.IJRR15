function ret = givens(n,d1,d2,t)
ret = eye(n);
ret([d1 d2],[d1 d2]) = [cos(t) sin(t); -sin(t) cos(t)];
end