function plot3dObstacles(obstacles, obsType)
% figure; hold on;

if obsType ~= 5
    for i=1:2:length(obstacles(:,1))
        a = obstacles(i,1);
        b = obstacles(i,2);
        c = obstacles(i,3);
        d = obstacles(i+1,1);
        e = obstacles(i+1,2);
        f = obstacles(i+1,3);
        
        plotcube([a-d b-e c-f],[d e f],0.5,[0.5 0.5 0.5])
    end
else
    [x,y,z] = sphere;
    for i=1:length(obstacles(:,1))
        p = obstacles(i,1:3);
        r = obstacles(i,4);
        surf(r*x+p(1),r*y+p(2),r*z+p(3),...
            'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeAlpha',0.5);
    end
end
end

