function [obstacles, optimal] = ...
    generateObstacles(obsLength, numObs, initial, goal, dim, obsType)

optimal = 0;
if obsType == 1
        rng(17);
    % randomly generate the obstacles lower corners
    corners = (1)*rand(numObs,dim)-obsLength/2;
    
    obstacles = zeros(numObs*2,dim);
    obstacles(1:2:numObs*2-1,:) = corners;
    obstacles(2:2:numObs*2,:) = corners + obsLength;
    
    % check to make sure the goal and initial states are not blocked
    validObs = sampleFree([initial; goal], obstacles, obsType);
    if isnan(sum(validObs(:,1)))
        obstacles = generateObstacles(obsLength, numObs, initial, goal, dim, obsType);
    end
    
elseif obsType == 2
    if dim == 2
        % given set
        obstacles = [0.434752 0.542681
            0.853504 0.947454
            0.543562 0.0691414
            0.945129 0.521948
            0.0529046 0.215034
            0.421858 0.617981];
    elseif dim == 3
        % 3d obstacle filled set
        obstacles = [0.434752 0.542681 -0.1
            0.853504 0.947454 1.1
            0.543562 0.0691414 -0.1
            0.945129 0.521948 1.1
            0.0529046 0.215034 -0.1
            0.421858 0.617981 1.1];
    end
elseif obsType == 3
    % maze
    obstacles = [0.7 0.3
        0.8 1.1 %
        0.3 0.3
        0.8 0.4 %
        0.3 0.4
        0.4 0.8 %
        -0.1 0.8
        0.25 0.85 %
        0.05 0.7
        0.3 0.75 %
        -0.1 0.6
        0.25 0.65 %
        0.05 0.5
        0.3 0.55 %
        -0.1 0.4
        0.25 0.45 %
        0.8 0.7
        0.95 0.75 %
        0.85 0.6
        1.1 0.65 %
        0.8 0.5
        0.95 0.55 %
        0.85 0.4
        1.1 0.45];
    optimal = 4.67112;
elseif obsType == 4
    % inside and outside obstacle
    obstacles = [
        -0.1 0.05
        0.9 0.07
        0.1 0.1
        1.1 0.12
        -0.1 0.15
        0.9 0.17
        0.1 0.2
        1.1 0.22
        -0.1 0.25
        0.9 0.27
        0.1 0.3
        1.1 0.32
        -0.1 0.55
        0.9 0.57
        0.1 0.6
        1.1 0.62
        -0.1 0.65
        0.9 0.67
        0.1 0.7
        1.1 0.72
        -0.1 0.75
        0.9 0.77
        0.1 0.8
        1.1 0.82
        0.1 0.85
        0.3 0.95
        0.4 0.85
        0.6 0.95
        0.7 0.85
        0.9 0.95];
elseif obsType == 5
    % randomly generate the obstacles center
    obstacles = [(1)*rand(numObs,dim) ones(numObs,1)*obsLength/2];
    if dim == 2
        obstacles = [0.7 0.7 0.2;...
            0.8 0.3 0.2;...
            0.29 0.5 0.2];
    elseif dim == 3
        obstacles = [0.7 0.7 0.7 0.3;...
            0.8 0.8 0.3 0.3;...
            0.3 0.2 0.5 0.3;...
            0.2 0.8 0.5 0.2];
    elseif dim == 10
        obstacles = [0.5*ones(1,dim) 0.45];
    end
    
    % check to make sure the goal and initial states are not blocked
    validObs = sampleFree([initial; goal], obstacles, obsType);
    if isnan(sum(validObs(:,1)))
        obstacles = generateObstacles(obsLength, numObs, initial, goal, dim, obsType);
    end
elseif obsType == 6 % bug trap
    obstacles = [0.25 0.33 % lower wall
        0.7 0.35
        0.25 0.63 % upper wall
        0.7 0.65
        0.68 0.33 % back wall
        0.7 0.65
        0.25 0.33 % front lower wall
        0.27 0.48
        0.25 0.52 % front upper wall
        0.27 0.65
        0.25 0.46 % mouth lower
        0.4 0.48
        0.25 0.52 % mouth upper
        0.4 0.54];
elseif obsType == 7 % recursive maze
    filename = ['maze/maze' num2str(dim)];
    obstacles = dlmread(filename);
elseif obsType == 8 % kinematic chain
    obstacles = [0.2 0.8
        0.4 1.0
        0.65 0.6
        0.75 1.0
        0.2 0.1
        0.3 0.4];
    if (isempty(sampleFree(initial,obstacles,obsType)) || ...
            isempty(sampleFree(goal,obstacles,obsType)))
        display('ERROR **Feas**: Initial or goal state is infeasible');
    end
elseif obsType == 9 % kinematic chain from Lavalle
        obstacles = [0.0 0.55
        0.15 0.6
        0.85 0.55
        1.0 0.6
        0.12 0.6
        0.15 0.8
        0.48 0.6
        0.52 0.8
        0.85 0.6
        0.88 0.8];
    if (isempty(sampleFree(initial,obstacles,obsType)) || ...
            isempty(sampleFree(goal,obstacles,obsType)))
        display('ERROR **Feas**: Initial or goal state is infeasible');
    end
end
end

