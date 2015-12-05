function [w,x,y,z]=euler2quat(rx,ry,rz,order)
% rx , ry & rz should be in radians
% rx, ry, rz are rotation around x,y & z axis
% order of rotations should be passed in 'order' parameter as string e.g. 'xyz' or
% 'xzy' etc.
% Quaterion is expressed as Q=[w x y z]
% w is scalar component of quaternion
% x,y,z are vector components of quaternion. For example q=w+xi+yj+zk



Qx=[cos(rx/2) sin(rx/2)  0          0];
Qy=[cos(ry/2) 0          sin(ry/2)  0];
Qz=[cos(rz/2) 0          0          sin(rz/2)];

if(strcmpi(order,'xyz')==1)  
%    disp('order is x-->y-->z');
    Q1=quaternionmul(Qy,Qz);  
    Q2=quaternionmul(Qx,Q1);  
    
elseif(strcmpi(order,'xzy')==1)  
 %   disp('order is x-->z-->y');
    Q1=quaternionmul(Qz,Qy);  
    Q2=quaternionmul(Qx,Q1);
    
elseif(strcmpi(order,'yxz')==1)
  %  disp('order is y-->x-->z');
    Q1=quaternionmul(Qx,Qz);  
    Q2=quaternionmul(Qy,Q1);
    
elseif(strcmpi(order,'yzx')==1)
   % disp('order is y-->z-->x');
    Q1=quaternionmul(Qz,Qx);  
    Q2=quaternionmul(Qy,Q1);
    
elseif(strcmpi(order,'zxy')==1)
    %disp('order is z-->x-->y');
    Q1=quaternionmul(Qx,Qy);  
    Q2=quaternionmul(Qz,Q1);
    
elseif(strcmpi(order,'zyx')==1) % same as matlab quatdemo
    %disp('order is z-->y-->x');
    Q1=quaternionmul(Qy,Qx);  
    Q2=quaternionmul(Qz,Q1); 
else
    %disp('could not recognized')    
end

    w=Q2(1,1);
    x=Q2(1,2);
    y=Q2(1,3);
    z=Q2(1,4);


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------    
% % Technical Reference: Ken Shoemake, "Animating Rotations with Quaternion Curves"
