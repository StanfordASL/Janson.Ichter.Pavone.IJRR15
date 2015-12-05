%% PLOT results
close all
% specify top level directory
filepath = '~/Desktop/ASL_Shared/Results/AA203/';
plotTitles = 'SE(3) Maze';
testcase = 'se3_';
successThreshold = 0.5;

% se2 sample counts
globalSampleCount = [9,64,100,144,196,256,324,576,729,900,1296,1764,2025,2304,2601,2916,3600,4356];
% se3 sample counts
globalSampleCount = [8,32,64,144,216,324,486,729,972,1296,1728,2304,3072,4096,5120,6400];
globalSampleCount = [9,64,100,144,196,256,324,576,729,900,1296,1764,2025,2304,2601,2916,3600,4356];

DataLabels = {'vanilla','parent'};
detSamples = [0 0]; % are the number of samples exact?
runCounts = [50 50];
legendLabels = {'FMT*','Parent-FMT*'};
styles = {'b-o','k-x'};
lineWidth = 1.4;

figCost = figure; hold on;
figTime = figure; hold on;
figSuccess = figure; hold on;
for dataIdx = 1:length(DataLabels)
    data = dlmread([filepath testcase DataLabels{dataIdx} '.txt'],',');
    runs = runCounts(dataIdx);
    counts = globalSampleCount;
    if(~detSamples(dataIdx))
        for i = 1:length(data(:,1)) % merge results
            bestFit = 0;
            miss = inf;
            for j = 1:length(counts)
                if abs(counts(j) - data(i,1)) < miss
                    miss = abs(counts(j) - data(i,1));
                    bestFit = counts(j);
                end
            end
            data(i,1) = bestFit;
        end
    end
    sampleCount = unique(data(:,1));
    sampleCount = globalSampleCount;
    cost = zeros(length(sampleCount),1);
    costStd = zeros(length(sampleCount),1);
    time = zeros(length(sampleCount),1);
    success = zeros(length(sampleCount),1);
    for i = 1:length(sampleCount)
        numSamples = sampleCount(i);
        success(i) = sum(numSamples == data(:,1));
        cost(i) = sum(data(numSamples == data(:,1),2))/success(i);
        costStd(i) = std(data(numSamples == data(:,1),2))/sqrt(success(i));
        time(i) = sum(data(numSamples == data(:,1),3))/success(i);
        success(i) = success(i)/runs;
    end
    cost(success < successThreshold) = nan;
    
    figure(figCost);
    errorbar(sampleCount,cost,costStd,styles{dataIdx},'LineWidth',lineWidth);
   
    figure(figTime);
    plot(time,cost,styles{dataIdx},'LineWidth',lineWidth);

    figure(figSuccess);
    plot(sampleCount,success,styles{dataIdx},'LineWidth',lineWidth);

end

% figure settings
figure(figCost);
axis([0 inf -inf inf]);
plotSettings(plotTitles,{'Sample Count','Solution Cost'},legendLabels,'NorthEast');
saveas(figCost,[filepath,'plots/',testcase,'cost.png']);

figure(figTime);
plotSettings(plotTitles,{'Time','Solution Cost'},legendLabels,'NorthEast');
saveas(figTime,[filepath,'plots/',testcase,'time.png']);

figure(figSuccess);
plotSettings(plotTitles,{'Sample Count','Success Rate'},legendLabels,'SouthEast');
saveas(figSuccess,[filepath,'plots/',testcase,'success.png']);