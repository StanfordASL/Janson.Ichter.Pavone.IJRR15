function [error] = checkErrors(dim, initial, goal, samplingType, planningType, obsType)
% function to check setup parameters and report any issues, exiting the
% program

error = 0;

% check that dimensions are all consistent
if dim ~= length(initial) || dim ~= length(goal)
    display('ERROR **Dims**: Dimensions of the initial point, goal point, and dim variable are not consistent!');
    error = 1;
end

% check that sampling technique and nearest neighbor functions are
% consistent
if (samplingType == 1 || samplingType == 4 || samplingType == 5) && (planningType > 1 && planningType ~= 7 && planningType ~= 8 && planningType ~= 9  && planningType ~= 10) && planningType ~= 6
    display('ERROR **Method**: Box nearest neighbor or lattice methods cannot be used with random sampling, change either samplingType or planningType');
    error = 1;
end

% check obstacle correctly configured
if dim ~= size(generateObstacles(0.1,1,initial,goal,dim,obsType),2) && obsType ~= 5 && obsType ~= 8 && obsType ~= 9
    display('ERROR **Obs**: Obstacle dimension is not consistent with dimension, check generateObstacles.m');
    error = 1;
end

% adaptive FMT not yet implemented for spherical obstacles or chains
if (obsType == 5 || obsType == 8 || obsType == 9) && planningType == 6
    error = 1;
    display('ERROR **Obs**: AdaptiveFMT does not yet support spherical obstacles');
end