clear all; close all;
dim = 10;
obsType = 9; % 1 for random, 2 for givenSet (only 2 or 3 dim), 3 for maze (2D),
% 4 for inside out, 5 for spherical, 6 for 2D bugtrap, 7 for recursive maze
coverage = 0.5;
smoothOn = 0;
saveFlag = 'v2';
resultsFolder = 'results_chain';
types = [
    1 8;
%     2 1;
    3 8;
    4 8;
%     5 1;
%     4 8;
%     6 8;
%     6 8;
    ];
% Each row is a [samplingType planningType]
% samplingType: 1 for rand (set planning type to 1), 2 for pert grid,
% 3 for grid, 4 for halton, 5 for hammersley
% planningType: 1 for normal NN, 2 for 4 move lattice, 3 for 8 move,
% 4 for 16 move 5 for index N, 6 for adaptive FMT, 7 for parent FMT, 8 for
% PRM*
successThreshold = 0.25;

plotTypes = {'Solution Cost vs Sample Count','Smooth Solution Cost vs Sample Count','Time vs Sample Count'...,
   'Solution Cost vs Time','Smooth Solution Cost vs Time','Success Rate vs Sample Count'};
lineColor = {'b',[0 0 0.8],'k',[0.9 0.7 0],'m','r'}; % will denote color of line based on samplingType
lineStyle = {'o-','+-.','x-.','*-.','s-','d--','^-','>-'}; % denote marker specifier based on planning Type

for plotType = 1:length(plotTypes)
    if(isempty(strfind(plotTypes{plotType},'Smooth')) || (smoothOn))
        h = figure; hold on;
        plotNum = 1;
        legendLables = cell(1,size(types,1)+1);
        for type = 1:size(types,1)
            % read data
            samp = types(type,1);
            plan = types(type,2);
            inputfile = genFilename(resultsFolder, dim, obsType, coverage,...
                samp, plan, saveFlag);
            results = dlmread(inputfile);
            splits = find(isnan(dlmread(inputfile)),2,'first');
            dataSplitEnd = splits(1)-1;
            stdSplitStart = splits(1)+1;
            stdSplitEnd = splits(2)-1;
            optSplit = splits(2) + 1;
            optimalCost = results(optSplit,1);
            
            success = results(1:dataSplitEnd,4); % success rate
            lowSamples = ones(size(success));
            lowSamples(success < successThreshold) = nan; % remove anything with less than 75% success
            
            % define colors and style
            color = lineColor{samp};
            style = lineStyle{plan};
            switch samp
                case 1 % random
                    color = 'b';
                    style = '-o';
                    switch plan
                        case 6
                            color = 'b';
                            style = '-x';
                        case 7
                            color = 'k';
                            style = '-x';
                    end
                case 2 % pert grid
                    color = 'k';
                    style = '-+';
                    switch plan
                        case 1
                            color = 'r';
                            style = '-+';
                        case 5
                            color = [0 0.6 0];
                            style = '-s';
                        case 6
                            color = 'r';
                            style = '-*';
                    end
                case 3 % grid
                    color = 'b';
                    style = '--s';
                    switch plan
                        case 1
                            color = 'k';
                            style = '-s';
                        case 5
                            color = [0 0.6 1];
                            style = '-x';
                        case 7
                            color = 'g';
                            style = '-x';
                        case 8
                            color = 'k';
                            style = '-s';
                        end
                case 4 % halton
                    color = 'g'; %[0 0.6 0];
                    style = '-.d';
                    switch plan
                        case 6
                            color = 'b';
                            style = '-x';
                        case 7
                            color = 'b';
                            style = '-x';
                    end
                case 5
                    color = [1 .5 0];
                    style = '-*';
                case 6 % hex grid
                    color = [1 .5 0];
                    color = 'k'
                    style = '-*';
                    style = '-+';
            end
            
            % create x, y, and error data
            axis([0 inf -inf inf]);
            switch plotType
                case 1 % cost vs samples
                    xdata = results(1:dataSplitEnd,1); % samples
                    ydata = results(1:dataSplitEnd,2).*lowSamples; % cost
                    edata = results(stdSplitStart:stdSplitEnd,2); % cost std
                    legendLoc = 'northeast';
                case 2 % smooth cost vs samples
                    xdata = results(1:dataSplitEnd,1); % samples
                    ydata = results(1:dataSplitEnd,3).*lowSamples; % smooth cost
                    edata = results(stdSplitStart:stdSplitEnd,3); % smooth std
                    legendLoc = 'best';
                case 3 % time vs samples
                    xdata = results(1:dataSplitEnd,1); % samples
                    ydata = results(1:dataSplitEnd,5).*lowSamples; % time
                    edata = results(stdSplitStart:stdSplitEnd,5); % time std
                    legendLoc = 'northwest';
                case 4 % cost vs time
                    xdata = results(1:dataSplitEnd,5); % time
                    ydata = results(1:dataSplitEnd,2).*lowSamples; % cost
                    edata = results(stdSplitStart:stdSplitEnd,2); % cost std
                    legendLoc = 'northeast';
                case 5 % smooth cost vs time
                    xdata = results(1:dataSplitEnd,5); % time
                    ydata = results(1:dataSplitEnd,3).*lowSamples; % smooth cost
                    edata = results(stdSplitStart:stdSplitEnd,3); % smooth std
                    legendLoc = 'best';
                case 6 % success rate vs samples
                    xdata = results(1:dataSplitEnd,1); % samples
                    ydata = results(1:dataSplitEnd,4); % success rate
                    edata = zeros(size(ydata));
                    legendLoc = 'best';
                    axis([0 max(xdata) 0 1.01]);
            end
            
            % plot
            lineWidth = 1.4;
            % either normally or with standard deviation bars
            errorbar(xdata,ydata,edata,style,'Color',color,'MarkerSize',8,'LineWidth',lineWidth);
%             plot(xdata,ydata,style,'Color',color,'MarkerSize',8,'LineWidth',lineWidth);
            % naming
            planningString = '';
            switch plan
                case 1
                    planningString = 'FMT*';
                case 2
                    planningString = '4 Move Lattice';
                case 3
                    planningString = '8 Move Lattice';
                case 4
                    planningString = '16 Move Lattice';
                case 5
                    planningString = 'Index  NN';
                case 6
                    planningString = 'Adaptive FMT*';
                case 7
                    planningString = 'Parent FMT*';
                case 8
                    planningString = 'PRM*';
            end
            samplingString = '';
            switch samp
                case 1
                    samplingString = 'i.i.d.';
                case 2
                    samplingString = 'Pert Grid';
                case 3
                    samplingString = 'Lattice';
                case 4
                    samplingString = 'Halton';
                case 5
                    samplingString = 'Random Lattice';
%                     samplingString = 'Hammersley';
                case 6
                    samplingString = 'Rand Rot Lattice';
            end
            % plot text
            legendLabels{plotNum} = [samplingString,', ',planningString];
            %         legendLabels{plotNum} = planningString;
            %         legendLabels{plotNum} = samplingString;
            plotNum = plotNum+1;
        end
        % plot optimal cost
        %         if(optimalCost~=0)
        %             plot(xlim(),[optimalCost optimalCost],'k-');
        %             legendLabels{plotNum} = 'Optimal';
        %         end
        % title and file naming
        switch obsType
            case 1
                obsDesc = num2str(round(coverage*100));
                titleString = [num2str(dim),'D',' Random, ',obsDesc,' Coverage'];
                saveString = [num2str(dim),'D',', Random, ',obsDesc,' Coverage, ',plotTypes{plotType}];
            case 2
                titleString = [num2str(dim),'D',' Rectangles'];
                saveString = [num2str(dim),'D',', Rectangles, ',plotTypes{plotType}];
            case 3
                titleString = [num2str(dim),'D',' Maze Obstacle'];
                saveString = [num2str(dim),'D',', Maze Obstacle, ',plotTypes{plotType}];
            case 4
                titleString = [num2str(dim),'D',' Vary Envi Obstacle'];
                saveString = [num2str(dim),'D',', Vary Envi Obstacle, ',plotTypes{plotType}];
            case 5
                obsDesc = num2str(round(coverage*100));
                titleString = [num2str(dim),'D',' Spherical'];
                saveString = [num2str(dim),'D',', Spherical, ',obsDesc,' Coverage, ',plotTypes{plotType}];
            case 6
                titleString = [num2str(dim),'D',' Bug Trap'];
                saveString = [num2str(dim),'D',', Bug Trap, ',plotTypes{plotType}];
            case 7
                titleString = [num2str(dim),'D',' Recursive Maze'];
                saveString = [num2str(dim),'D',', Recursive Maze, ',plotTypes{plotType}];
            case 8
                titleString = [num2str(dim),'D',' Kinematic Chain'];
                saveString = [num2str(dim),'D',', Kinematic Chain, ',plotTypes{plotType}];
            case 9
                titleString = [num2str(dim),'D',' Kinematic Chain Gap'];
                saveString = [num2str(dim),'D',', Kinematic Chain Gap, ',plotTypes{plotType}];
        end
        labelText = regexp(plotTypes{plotType},' vs ','split');
        plotSettings(titleString,{labelText{2},labelText{1}},legendLabels,legendLoc)
        % Save file
        saveString = strrep(saveString,',','_');
        saveString = strrep(saveString,' ','');
        saveas(h,[resultsFolder,'/plots/' saveString '.png']);
    end
end