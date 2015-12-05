%% PLOT path and tree
clear all; close all;

filepath = '~/Desktop/ASL_Shared/PianoMoversGrid/7D/48/FMT_1_Test/';
runNum = 14;
tree = dlmread([filepath 'tree_' num2str(runNum) '.txt']);
path = dlmread([filepath 'path_' num2str(runNum) '.txt']);

figure; 
dim1 = 1; dim2 = 2; dim3 = 3;
for i = 3:2:length(tree)
    plot3(tree(i:i+1,dim1),tree(i:i+1,dim2),tree(i:i+1,dim3),...
        'Color',0.4*ones(1,3));
    hold on;
end
plot3(path(2:length(path),dim1),path(2:length(path),dim2),path(2:length(path),dim3),'r');

%% PLOT results
close all
% specify top level directory
filepath = '~/Desktop/ASL_Shared/OMPL_Results/PianoMovers_Results/';
successThreshold = 0.75; 

% spatial grid results
SpatialGrid = dlmread([filepath 'all_costs_SpatialGrid.txt'],',');
runsSpatialGrid = 17;
samplesSpatialGrid = unique(SpatialGrid(:,1));
costSpatialGrid = zeros(length(samplesSpatialGrid),1);
costSpatialGridStd = zeros(length(samplesSpatialGrid),1);
timeSpatialGrid = zeros(length(samplesSpatialGrid),1);
successSpatialGrid = zeros(length(samplesSpatialGrid),1);
for i = 1:length(samplesSpatialGrid)
    numSamples = samplesSpatialGrid(i);
    successSpatialGrid(i) = sum(numSamples == SpatialGrid(:,1));
    costSpatialGrid(i) = sum(SpatialGrid(numSamples == SpatialGrid(:,1),2))/successSpatialGrid(i);
    costSpatialGridStd(i) = std(SpatialGrid(numSamples == SpatialGrid(:,1),2))/sqrt(successSpatialGrid(i));
    timeSpatialGrid(i) = sum(SpatialGrid(numSamples == SpatialGrid(:,1),3))/successSpatialGrid(i);
    successSpatialGrid(i) = successSpatialGrid(i)/runsSpatialGrid;
end
costSpatialGrid(successSpatialGrid < successThreshold) = nan;

RotGrid = dlmread([filepath 'all_costs_RotGrid.txt'],',');
runsRotGrid = 50;
sideCounts = [8 32 64 144 216 324 486 729 972 1296 1728 2304 3072 4096 5120 6400 8000];
for i = 1:length(RotGrid(:,1)) % merge results
    [~,idx] = min(abs(sideCounts-RotGrid(i,1)));
    RotGrid(i,1) = sideCounts(idx);
end
samplesRotGrid = unique(RotGrid(:,1));
costRotGrid = zeros(length(samplesRotGrid),1);
costRotGridStd = zeros(length(samplesRotGrid),1);
timeRotGrid = zeros(length(samplesRotGrid),1);
successRotGrid = zeros(length(samplesRotGrid),1);
for i = 1:length(samplesRotGrid)
    numSamples = samplesRotGrid(i);
    successRotGrid(i) = sum(numSamples == RotGrid(:,1));
    costRotGrid(i) = sum(RotGrid(numSamples == RotGrid(:,1),2))/successRotGrid(i);
    costRotGridStd(i) = std(RotGrid(numSamples == RotGrid(:,1),2))/sqrt(successRotGrid(i));
    timeRotGrid(i) = sum(RotGrid(numSamples == RotGrid(:,1),3))/successRotGrid(i);
    successRotGrid(i) = successRotGrid(i)/runsRotGrid;
end
costRotGrid(successRotGrid < successThreshold) = nan;

% random points results
RandomPoints = dlmread([filepath 'all_costs_RandomPoints.txt'],',');
runsRandomPoints = 50;
samplesRandomPoints = unique(RandomPoints(:,1));
costRandomPoints = zeros(length(samplesRandomPoints),1);
costRandomPointsStd = zeros(length(samplesRandomPoints),1);
timeRandomPoints = zeros(length(samplesRandomPoints),1);
successRandomPoints = zeros(length(samplesRandomPoints),1);
for i = 1:length(samplesRandomPoints)
    numSamples = samplesRandomPoints(i);
    successRandomPoints(i) = sum(numSamples == RandomPoints(:,1));
    costRandomPoints(i) = sum(RandomPoints(numSamples == RandomPoints(:,1),2))/successRandomPoints(i);
    costRandomPointsStd(i) = std(RandomPoints(numSamples == RandomPoints(:,1),2))/sqrt(successRandomPoints(i));
    timeRandomPoints(i) = sum(RandomPoints(numSamples == RandomPoints(:,1),3))/successRandomPoints(i);
    successRandomPoints(i) = successRandomPoints(i)/runsRandomPoints;
end
costRandomPoints(successRandomPoints < successThreshold) = nan;

% pert grid results
PertGrid = dlmread([filepath 'all_costs_PertGrid.txt'],',');
runsPertGrid = 50;
samplesPertGrid = unique(PertGrid(:,1));
costPertGrid = zeros(length(samplesPertGrid),1);
costPertGridStd = zeros(length(samplesPertGrid),1);
timePertGrid = zeros(length(samplesPertGrid),1);
successPertGrid = zeros(length(samplesPertGrid),1);
for i = 1:length(samplesPertGrid)
    numSamples = samplesPertGrid(i);
    successPertGrid(i) = sum(numSamples == PertGrid(:,1));
    costPertGrid(i) = sum(PertGrid(numSamples == PertGrid(:,1),2))/successPertGrid(i);
    costPertGridStd(i) = std(PertGrid(numSamples == PertGrid(:,1),2))/sqrt(successPertGrid(i));
    timePertGrid(i) = sum(PertGrid(numSamples == PertGrid(:,1),3))/successPertGrid(i);
    successPertGrid(i) = successPertGrid(i)/runsPertGrid;
end
costPertGrid(successPertGrid < successThreshold) = nan;

OffGrid = dlmread([filepath 'all_costs_OffGrid.txt'],',');
runsOffGrid = 20;
sideCounts = [8 32 64 144 216 324 486 729 972 1296 1728 2304 3072 4096 5120 6400 8000];
for i = 1:length(OffGrid(:,1)) % merge results
    [~,idx] = min(abs(sideCounts-OffGrid(i,1)));
    OffGrid(i,1) = sideCounts(idx);
end
samplesOffGrid = unique(OffGrid(:,1));
costOffGrid = zeros(length(samplesOffGrid),1);
costOffGridStd = zeros(length(samplesOffGrid),1);
timeOffGrid = zeros(length(samplesOffGrid),1);
successOffGrid = zeros(length(samplesOffGrid),1);
for i = 1:length(samplesOffGrid)
    numSamples = samplesOffGrid(i);
    successOffGrid(i) = sum(numSamples == OffGrid(:,1));
    costOffGrid(i) = sum(OffGrid(numSamples == OffGrid(:,1),2))/successOffGrid(i);
    costOffGridStd(i) = std(OffGrid(numSamples == OffGrid(:,1),2))/sqrt(successOffGrid(i));
    timeOffGrid(i) = sum(OffGrid(numSamples == OffGrid(:,1),3))/successOffGrid(i);
    successOffGrid(i) = successOffGrid(i)/runsOffGrid;
end
costOffGrid(successOffGrid < successThreshold) = nan;

% FullGrid = dlmread([filepath 'all_costs_FullGrid.txt'],',');

% plotting
styleSpatialGrid = 'g-x'; 
styleRotGrid = 'k-s';
styleRandomPoints = 'r-o';
stylePertGrid = 'b-+';
styleOffGrid = 'm->';
legendText = {'Spatial Grid','Grid','Random','Perturbed Grid','Offset Grid'};
fontSize = 16;
figure; hold on; set(gca,'FontSize',fontSize);
title('SE(3) Maze, Sample Count vs. Cost');
xlabel('Sample Count'); ylabel('Cost');
% plot(samplesSpatialGrid,costSpatialGrid,styleSpatialGrid);
% plot(samplesRotGrid,costRotGrid,styleRotGrid);
% plot(samplesRandomPoints,costRandomPoints,styleRandomPoints);
% plot(samplesPertGrid,costPertGrid,stylePertGrid);
% plot(samplesOffGrid,costOffGrid,styleOffGrid);
errorbar(samplesSpatialGrid,costSpatialGrid,costSpatialGridStd,styleSpatialGrid);
errorbar(samplesRotGrid,costRotGrid,costRotGridStd,styleRotGrid);
errorbar(samplesRandomPoints,costRandomPoints,costRandomPointsStd,styleRandomPoints);
errorbar(samplesPertGrid,costPertGrid,costPertGridStd,stylePertGrid);
errorbar(samplesOffGrid,costOffGrid,costOffGridStd,styleOffGrid);
legend(legendText,'Location','NorthEast');

figure; hold on; set(gca,'FontSize',fontSize);
title('SE(3) Maze, Time vs. Cost');
xlabel('Time'); ylabel('Cost');
plot(timeSpatialGrid,costSpatialGrid,styleSpatialGrid);
plot(timeRotGrid,costRotGrid,styleRotGrid);
plot(timeRandomPoints,costRandomPoints,styleRandomPoints);
plot(timePertGrid,costPertGrid,stylePertGrid);
plot(timeOffGrid,costOffGrid,styleOffGrid);
legend(legendText,'Location','NorthEast');

figure; hold on; set(gca,'FontSize',fontSize);
title('SE(3) Maze, SuccessRate vs. Cost');
xlabel('Sample Count'); ylabel('Success Rate');
plot(samplesSpatialGrid,successSpatialGrid,styleSpatialGrid);
plot(samplesRotGrid,successRotGrid,styleRotGrid);
plot(samplesRandomPoints,successRandomPoints,styleRandomPoints);
plot(samplesPertGrid,successPertGrid,stylePertGrid);
plot(samplesOffGrid,successOffGrid,styleOffGrid);
legend(legendText,'Location','SouthEast');
