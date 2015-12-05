% quaternion multiplication
function QMR=quaternionmul(Q0,Q1)
% Q0 and Q1 shoud be in this order
% Q0=[w0 x0 y0 z0] % w0 is scalar part, x0,y0,z0 are vector part
% Q1=[w1 x1 y1 z1] % w1 is scalar part, x1,y1,z1 are vector part
% Multiplication is not commutative in that the products Q0Q1 and Q1Q0 are
% not necessarily equal.
w0=Q0(1); x0=Q0(2); y0=Q0(3); z0=Q0(4); 
w1=Q1(1); x1=Q1(2); y1=Q1(3); z1=Q1(4); 

wr=(w0.*w1 - x0.*x1 - y0.*y1 - z0.*z1);
xr=(w0.*x1 + x0.*w1 + y0.*z1 - z0.*y1);
yr=(w0.*y1 - x0.*z1 + y0.*w1 + z0.*x1);
zr=(w0.*z1 + x0.*y1 - y0.*x1 + z0.*w1);

QMR=[wr xr yr zr]; % wr is scalar part, xr, yr, zr are vector part

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------

% Technical Reference: Ken Shoemake, "Animating Rotations with Quaternion Curves"
