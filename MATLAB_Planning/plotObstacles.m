function plotObstacles(obstacles, obsType)

if obsType ~= 5
    for i=1:2:length(obstacles(:,1))
        a = obstacles(i,1);
        b = obstacles(i,2);
        c = obstacles(i+1,1);
        d = obstacles(i+1,2);
        
        fill([a a c c], [d b b d], 'k','EdgeColor','k');
    end
else
    ang=0:0.01:2*pi;
    xp=cos(ang);
    yp=sin(ang);
    for i=1:length(obstacles(:,1))
        p = obstacles(i,1:2);
        r = obstacles(i,size(obstacles,2));
        fill(p(1)+r*xp,p(2)+r*yp,'k','EdgeColor','y');
    end
end
end

