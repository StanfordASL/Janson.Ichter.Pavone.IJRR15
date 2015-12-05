% Remove samples within obstacles

function [samples] = sampleFree(samples, obstacles, obsType)

if ~isempty(obstacles)
    dim = size(samples,2);
    
    notFree = false(length(samples),length(obstacles(:,1)));
    for i = 1:size(samples,1)
        if obsType == 5 % circular obstacles
            for j = 1:length(obstacles(:,1))
                r = obstacles(j,size(obstacles,2));
                p = obstacles(j,1:size(obstacles,2)-1);
                notFree(i,j) = 1;
                if norm(samples(i,:) - p) > r
                    notFree(i,j) = 0;
                end
            end
        elseif obsType == 8 || obsType == 9 % kinematic chain
            thetas = samples(i,:);
            linkLength = 1/dim/sqrt(2);
            theta = 0;
            w = ones(1,2)*0.5;
            for link = 1:dim
                theta = theta + thetas(link);
                v = w;
                w = v + [cos(theta) sin(theta)]*linkLength;
                if ~checkCollision(v,w,obstacles,1)
%                     plot([v(1) w(1)],[v(2) w(2)],'y','LineWidth',8); %
%                     plots collision
                    notFree(i,:) = 1;
                    break;
                end
            end
        else % rectangular obstacles
            for j = 1:2:length(obstacles(:,1))
                notFree(i,j) = 1;
                for k = 1:dim
                    notFree(i,j) = notFree(i,j) ...
                        && samples(i,k) > obstacles(j,k) ...
                        && samples(i,k) < obstacles(j+1,k);
                    if ~notFree(i,j)
                        break
                    end
                end
            end
        end
    end
    % mark invalid points with nan
    samples(sum(notFree,2)>0,:) = nan(sum(sum(notFree,2)>0),dim)-1;
end
% remove all invalid points
samples = samples(~isnan(samples(:,1)),:);