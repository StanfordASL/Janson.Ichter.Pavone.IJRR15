%% Collision Chcker
% Written by Brian Ichter
% Given a v and end state, check whether the resulting edge is free of
% collisions
%
% Given:
%   v- start state
%   w- end state
%   obstacles- obstacle set
%
% Return:
%   valid- true if the path between v and w is valid
%   freeDist- distance to obstacle path is invalid
%

function [valid, freeDist] = checkCollision(v, w, obstacles,obsType)
if obsType ~= 5 && obsType ~= 8 && obsType ~= 9 % obstacles are cubes
    freeDist = 0;
    
    bb_min = min(v,w);
    bb_max = max(v,w);
    
    for k = 1:2:length(obstacles(:,1))
        obs = obstacles(k:k+1,:);
        if ~BroadphaseValidQ(bb_min, bb_max,obs)
            if ~MotionValidQ(v, w, obs)
                valid = 0;
                return;
            end
        end
    end
elseif obsType == 5 % obstacles are spheres
    if ~all(v == w)
        for k = 1:size(obstacles,1)
            r = obstacles(k,size(obstacles,2));
            p = obstacles(k,1:size(obstacles,2)-1);
            
            v2w = w-v;
            p2v = v-p;
            recipNorm = 1/norm(v2w);
            
            t = -dot(p2v, v2w)*recipNorm^2;
            dist = norm((p2v)+t*v2w);
            
            if(dist < r && t > 0.0 && t < 1.0)
                valid = 0;
                return;
            end
        end
    end
else % kinematic chain
    dim = size(w,2);
    linkLength = 1/dim/sqrt(2);
    
    % determine number of collision checks
    disc = 0.05;
    init = ones(1,2)*0.5;
    th_joint_init = 0;
    final = ones(1,2)*0.5;
    th_joint_final = 0;
    for link = 1:dim
        th_joint_init = th_joint_init + v(link);
        init = init + linkLength*[cos(th_joint_init) sin(th_joint_init)];
        th_joint_final = th_joint_final + w(link);
        final = final + linkLength*[cos(th_joint_final) sin(th_joint_final)];
    end
    dist = norm(init-final);
    itrs = ceil(dist/disc);
    du = mod(mod(w-v,2*pi)+pi,2*pi)-pi;
    
%     figure; hold on;
%     plotChain(v,'b');
%     plotChain(w,'g');
%     plotObstacles(obstacles,obsType);
    for itr = 1:itrs-1 % iterate through discretized configurations
        u = v + du*itr/itrs;
        v_joint = ones(1,2)*0.5; % initial location for all chains is center
        th_joint = u(1);
        w_joint = v_joint + [cos(th_joint) sin(th_joint)]*linkLength;
%         plotChain(u,'y');
        if ~checkCollision(v_joint, w_joint, obstacles, 1)
%             plotChain(u,'r');
%             plot([v_joint(1) w_joint(1)],[v_joint(2) w_joint(2)],'c','LineWidth',6);
            valid = 0;
            return;
        end
        for link = 2:dim % iterate through each link of the chain
            v_joint = w_joint;
            th_joint = th_joint + u(link);
            w_joint = v_joint + [cos(th_joint) sin(th_joint)]*linkLength;
            % call the above function (for two points) to see if a link is free
            if ~checkCollision(v_joint, w_joint, obstacles, 1)
%                 plotChain(u,'r');
%                 plot([v_joint(1) w_joint(1)],[v_joint(2) w_joint(2)],'c','LineWidth',6);
                valid = 0;
                return;
            end
        end
    end
end
valid = 1;

    function [isValid] = BroadphaseValidQ(bb_min, bb_max, obs)
        for i = 1:size(obs,2)
            if bb_max(i) <= obs(1,i)  || obs(2,i) <= bb_min(i)
                isValid = 1;
                return
            end
        end
        isValid = 0;
    end

% validMat = or(bb_max < obs(1,:),obs(2,:) < bb_min);
% if sum(validMat) > 0
%     isValid = 1;

    function [isValid] = MotionValidQ(v, w, obs)
        v_to_w = w-v;
        corner = v < obs(1,:);
        lambdas = (corner.*obs(1,:) + ~corner.*obs(2,:) - v)./v_to_w+5*eps;
        for i = 1:size(obs,2)
            if FaceContainsProjection(v, v_to_w, lambdas(i), i, obs)
                isValid = 0;
                minLambda = min(abs(lambdas));
                freeDist = norm(v_to_w)*minLambda;
                %                 figure; hold on;
                %                 plot([v(1) w(1)],[v(2) w(2)]);
                %                 plot(v(1),v(2),'g*');
                %                 plot(v(1)+v_to_w(1)*minLambda,v(2)+v_to_w(2)*minLambda,'rs');
                %                 plotObstacles(obs);
                return
            end
        end
        isValid = 1;
    end

    function [isValid] = FaceContainsProjection(v, v_to_w, lambda, j, obs)
        for i = 1:size(obs,2)
            if i ~= j && ~(obs(1,i) <= v(i) + v_to_w(i)*lambda && ...
                    v(i) + v_to_w(i)*lambda <= obs(2,i))
                isValid = 0;
                return
            end
        end
        isValid = 1;
    end
end