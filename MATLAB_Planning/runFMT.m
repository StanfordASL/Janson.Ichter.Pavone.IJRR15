%% FMT Simulation
% Written by Brian Ichter
% Main running file for FMT simulations
%% Simulation
% cd('/Users/brianichter/Desktop/ASL_Shared/Code/MATLAB_Planning');

clear all; close all;
format compact;
disp(['____________ New Run Started at ',datestr(now),' _______________']);

% Setup Params
global plotOn samplingType planningType indexMat2Array array2IndexMat gridCount dim;

% Run parameters for simulation
% dim) dimension
% sampling) 1 for rand (set planning type to 1), 2 for pert grid, 3 for
% grid, 4 for halton, 5 for hammrsley, 6 for hex 2D grid
% planning) 1 for normal NN, 2 for 4 move lattice, 3 for 8 move, 4 for 16 move
% 5 for index N, 6 for adaptive FMT, 7 for parent-FMT, 8 for PRM, 9 for
% bucket FMT, 10 for gpuFMT
% runCount) number of runs
% obsType) 1 for random, 2 for givenSet (only 2 or 3 dim), 3 for maze (2D),
% 4 for inside out, 5 for spherical, 6 for 2D bugtrap, 7 for recursive maze
% 8 for kinematic chain, 9 for kinematic chain gap

parameters = [...
    % dim     % sampling      % planning     % runCount     % obsType   % dummy
%     2         1               9              1               2;
%     2         1               1              1               2;
    3         1               1              1               2;
%     10         4               8              1               9;
%     10         1               8              50               9;
%     10         2               8              10               9;
%     2         3               8              1               3;
%     8         3               1              1               9;
    ]

initialLoc = 0.05; % position on the diagonal from 0 to 1
goalLoc = 0.95; % position on the diagonal from 0 to 1

sampleCounts = []; % set in for loop depending on dim;

% obstacle settings for type 1, random
obsLength = 0.2; % length each obstacle side (will be a n-dimensional shape with equal side lengths) or diameter
coverage = 0.6; % percentage of obstacle coverage (only for random obstacles)

% toggles and saving, 0 for off, 1 for on
saveFlag = ''; % empty ('') to not save
resultsFolder = 'results_bucket';
plotOn = 0;
smoothOn = 0;
forExcel = 0; % output individual columns for fast pasting to excel
rng('default'); % setup seed
single = true; % single data type, use for GPU computing

for param = 1:size(parameters,1)
    dim = parameters(param,1);
    if dim == 2
        sampleCounts = [2 3 5 7 8 12 15 18 21 24 28 32 36 40 45 50 55 60 70 80].^dim-1;
        sampleCounts = [1500];
    elseif dim == 3
        sampleCounts = [2 3 4 5 6 7 8 9 10 11 12 13 14 15].^dim-1;
        sampleCounts = [6 8 10 14].^dim-1;
        sampleCounts = 1500;
    elseif dim == 4
        sampleCounts = [2 3 4 4 5 5 6 6 7 7 8;
            2 3 4 4 5 5 6 6 7 7 8;
            2 3 3 4 4 5 5 6 6 7 7;
            2 3 3 4 4 5 5 6 6 7 7];
        sampleCounts = prod(sampleCounts);
%         sampleCounts = 100;
    elseif dim == 5
        sampleCounts = [2 3 3 4 4 4 4 5 5 5 5 5;
            2 3 3 3 4 4 4 4 5 5 5 5;
            2 2 3 3 4 4 4 4 4 5 5 5;
            2 2 3 3 3 4 4 4 4 4 5 5;
            2 2 2 3 3 3 4 4 4 4 4 5];
    elseif dim == 6
        sampleCounts = [2 3 3 3 4 4 4 4 5;
            2 3 3 3 4 4 4 4 4;
            2 2 3 3 3 4 4 4 4;
            2 2 3 3 3 4 4 4 4;
            2 2 2 3 3 3 4 4 4;
            2 2 2 3 3 3 3 4 4];
        sampleCounts = [2 2 2 2 2 2]';
    elseif dim == 7
        sampleCounts = [2 2 2 3 3 3 3 3 3 3 4 4;
            2 2 2 2 3 3 3 3 3 3 3 4;
            2 2 2 2 2 3 3 3 3 3 3 3;
            2 2 2 2 2 2 3 3 3 3 3 3;
            2 2 2 2 2 2 2 3 3 3 3 3;
            1 2 2 2 2 2 2 2 3 3 3 3;
            1 1 2 2 2 2 2 2 2 3 3 3];
    elseif dim == 8
        sampleCounts = [2 2 2 2 2 2 2 1;
            2 2 2 2 2 2 2 2;
            3 3 2 2 2 2 2 2;
            3 3 3 3 2 2 2 2;
            3 3 3 3 3 2 2 2;
            3 3 3 3 3 3 2 2;
            3 3 3 3 3 3 3 2;
            3 3 3 3 3 3 3 3;
            4 3 3 3 3 3 3 3;
            4 4 3 3 3 3 3 3]';
%             4 4 4 3 3 3 3 3;
%             4 4 4 4 3 3 3 3]';
%         sampleCounts = [4 4 4 4 4 3 3 3;
%             4 4 4 4 4 4 3 3;
%             ]';
        sampleCounts  = [3 3 3 3 2 2 2 2]';
    elseif dim == 10
        sampleCounts = [2 2 2 2 2 1 1 1 1 1
            2 2 2 2 2 2 1 1 1 1;
            2 2 2 2 2 2 2 1 1 1;
            2 2 2 2 2 2 2 2 1 1;
            2 2 2 2 2 2 2 2 2 1;
            2 2 2 2 2 2 2 2 2 2;
            3 2 2 2 2 2 2 2 2 2;
            3 3 2 2 2 2 2 2 2 2;
            3 3 3 2 2 2 2 2 2 2;
            3 3 3 3 2 2 2 2 2 2;
%             3 3 3 3 3 2 2 2 2 2;
%             3 3 3 3 2 2 2 2 2 2;
            ]';
        sampleCounts = [3 3 3 3 2 2 2 2 2 2]';
    elseif dim == 12
        sampleCounts = [2 2 2 2 2 2 1 1 1 1 1 1;
            2 2 2 2 2 2 2 1 1 1 1 1;
            2 2 2 2 2 2 2 2 1 1 1 1;
            2 2 2 2 2 2 2 2 2 1 1 1;
            2 2 2 2 2 2 2 2 2 2 1 1;
            2 2 2 2 2 2 2 2 2 2 2 1;
            2 2 2 2 2 2 2 2 2 2 2 2;
            3 2 2 2 2 2 2 2 2 2 2 2;
            3 3 2 2 2 2 2 2 2 2 2 2;
            3 3 3 2 2 2 2 2 2 2 2 2;
            ]';
        sampleCounts = [
            2 2 2 2 2 2 2 1 1 1 1 1;
            2 2 2 2 2 2 2 2 2 1 1 1]';
    else
        sampleCounts = [30 500 1000 3000];
    end
    samplingType = parameters(param,2);
    planningType = parameters(param,3);
    runCount = parameters(param,4);
    obsType = parameters(param,5);
    disp(['_________ Starting: dim = ',num2str(dim),', samp = ',num2str(samplingType),...
        ', plan = ',num2str(planningType),' , runs = ', num2str(runCount),' , obsType = ', num2str(obsType),' _________'])
    
    initial = initialLoc*ones(1,dim);
    goal = goalLoc*ones(1,dim); % only supports ending at [1.0 ... 1.0];
    
    % clear gpu if gpuFMT
    if planningType == 10
        d = gpuDevice;
    end
    
    if obsType == 8 % kinematic chain starting location
        % normal
        initial = zeros(1,dim);
        %         goal = zeros(1,dim);
        kink = ceil(dim/2)+1;
        initial(kink) = pi/2;
        %         goal(1) = -pi;
        %         goal(kink) = -pi/2;
        goal = -ones(1,dim)*2*pi/(dim+1);
        goal(1) = pi;
    elseif obsType == 9
        % Lavalle problem
        initial = zeros(1,dim);
        initial(1) = pi;
        initial(3) = -pi/2;
        goal = zeros(1,dim);
        goal(3) = pi/2;
    end
    
    if coverage > 0
        if obsType ~= 5
            obsVol = obsLength^dim;
        else
            obsVol = pi^(dim/2)*(obsLength/2)^dim/gamma(dim/2+1);
        end
        numObs = floor(coverage/obsVol);
    end
    
    error = checkErrors(dim, initial, goal, samplingType, planningType, obsType);
    if error
        return
    end
    
    % create arrays for storing data
    runCosts = zeros(size(sampleCounts,2),1);
    runCostsSmooth = zeros(size(sampleCounts,2),1);
    runTimes = zeros(size(sampleCounts,2),1);
    runCostsStd = zeros(size(sampleCounts,2),1);
    runCostsSmoothStd = zeros(size(sampleCounts,2),1);
    runTimesStd = zeros(size(sampleCounts,2),1);
    runCostsAll = zeros(size(sampleCounts,2),runCount);
    runCostsSmoothAll = zeros(size(sampleCounts,2),runCount);
    runTimesAll = zeros(size(sampleCounts,2),runCount);
    successRate = zeros(size(sampleCounts,2),1);
    
    % iterate through various sampleCounts
    for count = 1:size(sampleCounts,2)
        currentSampleCount = sampleCounts(1,count);
        sides= zeros(1,dim);
        if(size(sampleCounts,1) > 1)
            sides = sampleCounts(:,count);
            currentSampleCount = prod(sides);
        end
        disp(currentSampleCount);
        successes = 0;
        for runs = 1:runCount
            rng(runs);
            clear nodes; % required for timing
            startTime = now;
            
            if(size(sampleCounts,1) == 1)
                num = floor(sampleCounts(count).^(1/dim))+1;
                sampleCount = num^dim;
            else
                num = 0;
                sampleCount = currentSampleCount;
            end
            % Generate samples
            if(samplingType == 1)
                samples = rand(sampleCount,dim);
            elseif(samplingType == 2)
                if(size(sampleCounts,1) == 1)
                    samples = generatePertGridSamples(dim,num);
                else
                    samples = generateWeightedPertGridSamples(sides');
                end
            elseif(samplingType == 3)
                %                 samples = generateGridSamples(dim,num);
                if(size(sampleCounts,1) == 1)
                    samples = generateRotatedGridSamples(dim,num);
%                     samples = generateGridSamples(dim,num);
                else
                    if obsType == 8 || obsType == 9
                        samples = generateChainWeightedRotatedGridSamples(sides');
%                         samples = generateWeightedGridSamples(sides');
                    else
                        samples = generateWeightedRotatedGridSamples(sides');
                        %                     samples = generateWeightedGridSamples(sides');
                    end
                end
            elseif(samplingType == 4)
                samples = generateHaltonSamples(dim,sampleCount);
            elseif(samplingType == 5)
                samples = generateHammersleySamples(dim,sampleCount);
            elseif(samplingType == 6)
                samples = generateRotatedHexGridSamples(dim,sampleCount);
                %                 samples = generateHexGridSamples(dim,sampleCount);
            end
            
            if obsType == 8 || obsType == 9
                samples = samples*2*pi-pi;
            end
            
            %             disp(['Sample count pre sample free = ' num2str(size(samples,1))])
            
            if samplingType > 1
                [~,goalIdx] = min(sum(((samples - repmat(goal,size(samples,1),1))).^2,2));
                samples(goalIdx,:) = goal; % add the goal state to the samples
            else
                goalIdx = size(samples,1);
                samples(goalIdx,:) = goal; % add the goal state to the samples
            end
            
            % create translation matrix from array to matrix
            if planningType ~= 1 && planningType ~= 6 && planningType ~= 8 ...
                    && planningType ~= 7 && planningType ~= 9 && planningType ~= 10
                gridCount = round((size(samples,1))^(1/dim));
                matSize = cell(1,dim);
                for i = 1:dim
                    matSize{i} = gridCount;
                end
                indexMat2Array = reshape(find(samples(:,1)+1),matSize{:});
                subs = cell(1,dim);
                [subs{:}] = ind2sub(size(indexMat2Array),find(samples(:,1)+1));
                array2IndexMat = [subs{:}]';
            end
            
            
            % generate random obstacles
            [obstacles, optimal] = generateObstacles(obsLength, numObs, initial, goal, dim, obsType);
            
            % remove samples that are within the obstacle
            samples = sampleFree(samples,obstacles, obsType);
            
            %             disp(['Sample count post sample free = ' num2str(size(samples,1))])
            sampleCount = size(samples,1);
            % set ball radius for search
            eta = 0.1; %sqrt(dim)-0.9; %0.1;
            gammaval = (1+eta)*2*(1/dim)^(1/dim)*((1-coverage)/calculateUnitBallVolume(dim))^(1/dim);
            r = gammaval*(log(sampleCount)/sampleCount)^(1/dim);
            if planningType == 2
                r = 0; % r has no effect
            elseif planningType == 3
                r = 1/(sampleCount^(1/dim)); % grab nearest 8 nodes
            elseif planningType == 4
                r = 2/(sampleCount^(1/dim)); % grab nearest 20 nodes
            elseif planningType == 7
                %                 r = 1.5/(sampleCount^(1/dim));
            elseif planningType == 8 % PRM radius
                gammaval = (1+eta)*2*(1+1/dim)^(1/dim)*((1-coverage)/calculateUnitBallVolume(dim))^(1/dim);
                r = gammaval*(log(sampleCount)/sampleCount)^(1/dim);
                
                %                 if parameters(param,6) == 2
                %                     r = 8*gammaval*(sampleCount^(-0.5)/sampleCount)^(1/dim);
                %                     saveFlag = 'nmHalf_high';
                %                 end
                %                 if parameters(param,6) == 1
                %                     r = gammaval*(log(log(sampleCount))/sampleCount)^(1/dim);
                %                     saveFlag = 'loglogn_high';
                %                 end
            end
            
                
            if obsType == 8 || obsType == 9
                r = 1.15*sqrt(dim)*pi*r;
                if planningType == 8
                    r = r*0.75;
                end
            end
            
            % define potential function
%             potential = ones(size(samples));
            
            if single
                type = 'single';
                samples = cast(samples,type);
                initial = cast(initial,type);
                goal = cast(goal,type);
                obstacles = cast(obstacles,type);
                r = cast(r,type);
            end
            
            % Run FMT
            if planningType == 6 % adaptive FMT
                [path, nodes, costAll] = ...
                    adaptiveFMT(samples, initial, goal, obstacles, r);
            elseif planningType == 8 % PRM
                [path, nodes, costAll] = ...
                    PRM(samples, initial, goal, obstacles, r, obsType);
            elseif planningType == 9 % bucket FMT
                [path, nodes, costAll] = ...
                    z_bucketFMT(samples, initial, goal, obstacles, r, obsType);
            elseif planningType == 10 % gpuFMT
                [path, nodes, costAll] = ...
                    gpuFMT(samples, initial, goal, obstacles, r, obsType);
            else % FMT
                [path, nodes, costAll] = ...
                    FMT(samples, initial, goal, obstacles, r, obsType);
            end
            
            if plotOn % plot total tree
                % create figure
                figure; hold on;
                maxCost = max(costAll);
%                 maxCost = 40;
                if obsType == 8 || obsType == 9
                    axis([0.5-1/sqrt(2) 0.5+1/sqrt(2) 0.5-1/sqrt(2) 0.5+1/sqrt(2)]);
                    plotObstacles(obstacles, obsType)
                    
                    for i = 1:size(samples,1)
                        color = [nodes(i).cost/maxCost 0 1-nodes(i).cost/maxCost];
                        lineWidth = [];
                        if(isnan(nodes(i).cost));
                            color = [0 1 0];
                            lineWidth = 1;
                        end
                        plotChain(nodes(i).loc,color,lineWidth);
                    end
                    
                    plotChain(initial,'k',[],'-.');
                    plotChain(goal,'k',[],'-.');
                else
                    axis([0 1 0 1]);
                    %                 axis([min(samples(:,1)) max(samples(:,1))...
                    %                     min(samples(:,2)) max(samples(:,2))]);
                    
                    % plot obstacles
                    plotObstacles(obstacles, obsType)
                    
                    % plot tree
                    for i = 1:size(samples,1)
                        color = [nodes(i).cost/maxCost 0 1-nodes(i).cost/maxCost];
                        if(isnan(nodes(i).cost));
                            color = [0 1 0];
                        end
                        if(~isempty(nodes(i).Prev))
                            plot([nodes(i).loc(1) nodes(i).Prev.loc(1)],...
                                [nodes(i).loc(2) nodes(i).Prev.loc(2)],'-','Color',0.3*ones(1,3));
                        end
                        plot(nodes(i).loc(1),nodes(i).loc(2),...
                            '.','MarkerEdgeColor',color,'MarkerSize',15)
                    end
                end
            end
            
            % update goalIdx
            if samplingType > 1
                [~,goalIdx] = min(sum(((samples - repmat(goal,size(samples,1),1))).^2,2));
            else
                goalIdx = size(samples,1);
            end
            % if successful, plot solution and update results
            if(~isnan(costAll(goalIdx)))
                successes = successes + 1;
                solutionPath = nodes(goalIdx);
                if obsType ~= 8 && obsType ~= 9
                    updateCost(solutionPath);
                end
                runCostsAll(count,runs) = solutionPath.cost;
                runTimesAll(count,runs) = (now - startTime)*10^5;
                
                solutionPath = nodes(goalIdx);
                if plotOn % plot original solutions
                    if obsType == 8 || obsType == 9
                        figure; hold on; axis([0.5-1/sqrt(2) 0.5+1/sqrt(2) 0.5-1/sqrt(2) 0.5+1/sqrt(2)]);
                        plotObstacles(obstacles,obsType);
                        plotNode = solutionPath;
                        solnCost = plotNode.cost;
                        while(~isempty(plotNode) && ~isempty(plotNode.Prev))
                            color = [plotNode.cost/solnCost 1-plotNode.cost/solnCost 0];
                            plotChain(plotNode.loc,color);
                            plotNode = plotNode.Prev;
                        end
                        color = [plotNode.cost/solnCost 1-plotNode.cost/solnCost 0];
                        plotChain(plotNode.loc,color);
                    else
                        plotNode = solutionPath;
                        while(~isempty(plotNode) && ~isempty(plotNode.Prev))
                            plot([plotNode.loc(1) plotNode.Prev.loc(1)],...
                                [plotNode.loc(2) plotNode.Prev.loc(2)],'g-','LineWidth',3);
                            plotNode = plotNode.Prev;
                        end
                    end
                end
                
                if smoothOn
                    % smooth solution
                    solutionPath = smoothPath(solutionPath, obstacles, obsType);
                    updateCost(solutionPath);
                    runCostsSmoothAll(count,runs) = solutionPath.cost;
                    if plotOn % plot smoothed solution
                        if obsType == 8 || obsType == 9
                            
                        else
                            plotNode = solutionPath;
                            while(~isempty(plotNode) && ~isempty(plotNode.Prev))
                                plot([plotNode.loc(1) plotNode.Prev.loc(1)],...
                                    [plotNode.loc(2) plotNode.Prev.loc(2)],'m-','LineWidth',3);
                                plotNode = plotNode.Prev;
                            end
                        end
                    end
                end
                
            end
            
            if plotOn
                outputfile = strrep(genFilename([resultsFolder,'/plots'], dim, obsType, coverage,...
                    samplingType, planningType, saveFlag),'csv','png');
                saveas(gcf,outputfile);
            end
            
            % plot 3d results
            if plotOn && dim == 3
                figure;
                if obsType == 8 || obsType == 9
                    
                else
                    plot3(initial(1),initial(2),initial(3),'g*');
                    axis([min(samples(:,1)) max(samples(:,1))...
                        min(samples(:,2)) max(samples(:,2))...
                        min(samples(:,3)) max(samples(:,3))]);
                    hold on;
                    
                    % plot obstacles
                    plot3dObstacles(obstacles, obsType)
                    
                    % plot tree
                    for i = 1:size(samples,1)
                        color = [nodes(i).cost/max(costAll) 0 1-nodes(i).cost/max(costAll)];
                        if(isnan(nodes(i).cost));
                            color = [0 1 0];
                        end
                        color = round(color,5); % protect against floating point errors
                        if(~isempty(nodes(i).Prev))
                            plot3([nodes(i).loc(1) nodes(i).Prev.loc(1)],...
                                [nodes(i).loc(2) nodes(i).Prev.loc(2)],...
                                [nodes(i).loc(3) nodes(i).Prev.loc(3)],'r-');
                        end
                        plot3(nodes(i).loc(1),nodes(i).loc(2),nodes(i).loc(3),...
                            '.','MarkerEdgeColor',color)
                    end
                    
                    if(~isnan(costAll(goalIdx)))
                        solutionPath = nodes(goalIdx);
                        plotNode = solutionPath;
                        while(~isempty(plotNode) && ~isempty(plotNode.Prev))
                            plot3([plotNode.loc(1) plotNode.Prev.loc(1)],...
                                [plotNode.loc(2) plotNode.Prev.loc(2)],...
                                [plotNode.loc(3) plotNode.Prev.loc(3)],'m-','LineWidth',3);
                            plotNode = plotNode.Prev;
                        end
                    end
                end
            end
            
            % plot figure for adaptive FMT*
            if plotOn && planningType == 6
                % create figure
                figure; hold on;
                if obsType == 8 || obsType == 9
                    
                else
                    title('R plot');
                    axis([min(samples(:,1)) max(samples(:,1))...
                        min(samples(:,2)) max(samples(:,2))]);
                    
                    % plot obstacles
                    plotObstacles(obstacles, obsType)
                    
                    % plot tree
                    maxR = 0;
                    minR = inf;
                    for i = 1:size(samples,1)
                        if nodes(i).r > maxR
                            maxR = nodes(i).r;
                        end
                        if nodes(i).r < minR
                            minR = nodes(i).r;
                        end
                    end
                    
                    for i = 1:size(samples,1)
                        color = [(nodes(i).r-minR)/(maxR-minR) ...
                            0 ...
                            1-(nodes(i).r-minR)/(maxR-minR)];
                        if(isnan(nodes(i).cost));
                            color = [0 1 0];
                        end
                        if(~isempty(nodes(i).Prev))
                            plot([nodes(i).loc(1) nodes(i).Prev.loc(1)],...
                                [nodes(i).loc(2) nodes(i).Prev.loc(2)],'-','Color',0.3*ones(1,3));
                        end
                        plot(nodes(i).loc(1),nodes(i).loc(2),'.','MarkerEdgeColor',color,'MarkerSize',18);
                    end
                    
                    if(~isnan(costAll(goalIdx)))
                        plotNode = solutionPath;
                        while(~isempty(plotNode) && ~isempty(plotNode.Prev))
                            plot([plotNode.loc(1) plotNode.Prev.loc(1)],...
                                [plotNode.loc(2) plotNode.Prev.loc(2)],'m-','LineWidth',3);
                            plotNode = plotNode.Prev;
                        end
                    end
                    
                    outputfile = strrep(genFilename([resultsFolder,'/plots'], dim, obsType, coverage,...
                        RandomPoints, planningType, [num2str(saveFlag),'Rplot']),'csv','png');
                    saveas(gcf,outputfile);
                end
            end
            disp([num2str(runs),...
                ' -- time = ',num2str(runTimesAll(count,runs)),...
                ' -- cost = ',num2str(runCostsAll(count,runs))]);
        end
        runCosts(count) = sum(runCostsAll(count,:))/successes;
        runCostsSmooth(count) = sum(runCostsSmoothAll(count,:))/successes;
        runTimes(count) = sum(runTimesAll(count,:))/successes; % in seconds
        successRate(count) = successes/runCount;
        
        runCostsStd(count) = std(runCostsAll(count,runCostsAll(count,:) ~= 0))/sqrt(successes);
        runCostsSmoothStd(count) = std(runCostsSmoothAll(count,runCostsAll(count,:) ~= 0))/sqrt(successes);
        runTimesStd(count) = std(runTimesAll(count,runCostsAll(count,:) ~= 0))/sqrt(successes);
    end
    
    % Update sample counts to real sample counts
    if size(sampleCounts,1) == 1
        printSampleCounts = (floor(sampleCounts'.^(1/dim))+1).^dim;
    else
        printSampleCounts = prod(sampleCounts)';
    end
    
    % Print results
    disp(['______________ New Run Results at ',datestr(now),' _____________']);
    
    format short g
    results = [printSampleCounts runCosts runCostsSmooth successRate runTimes];
    X = gallery('uniformdata',size(results),0);
    disp('  Num Samples     Run Cost  Smooth Cost      Success     Run Time');
    disp(results)
    disp('Standard Errors');
    resultsStd = [printSampleCounts runCostsStd runCostsSmoothStd successRate runTimesStd];
    disp(resultsStd);
    
    disp(['Optimal cost = ',num2str(optimal)]);
    disp(['Run conditions: runCount = ',num2str(runCount),...
        ' , dim = ',num2str(dim),...
        ' , samplingType = ',num2str(samplingType),...
        ' , planningType = ',num2str(planningType)]);
    disp(['num obstacles = ',num2str(numObs),...
        ' , obstacle length = ',num2str(obsLength),...
        ' , obstacle type = ',num2str(obsType)]);
    
    if ~plotOn && ~isempty(saveFlag)
        outputfile = genFilename(resultsFolder, dim, obsType, coverage,...
            samplingType, planningType, saveFlag);
        dlmwrite(outputfile,results,'delimiter',',');
        dlmwrite(outputfile,nan(1,size(results,2)),'delimiter',',','-append');
        dlmwrite(outputfile,resultsStd,'delimiter',',','-append');
        dlmwrite(outputfile,nan(1,size(results,2)),'delimiter',',','-append');
        dlmwrite(outputfile,optimal,'delimiter',',','-append');
        dlmwrite(outputfile,[runCount dim ...
            samplingType planningType],'delimiter',',','-append');
        dlmwrite(outputfile,[numObs obsLength obsType]...
            ,'delimiter',',','-append');
    end
    
    if forExcel
        printSampleCounts
        runCosts
        runCostsSmooth
        successRate
        runTimes
        runCostsStd
        runCostsSmoothStd
        runTimesStd
    end
    format short
    
    % clear gpu if gpuFMT
    if planningType == 10
        d = gpuDevice;
        reset(d)
    end
end