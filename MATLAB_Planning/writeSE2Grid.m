%% To generate rotated grid in OMPL for the SE(2) Maze
% Uses rotated pert grids for the spatial and a random angle
clear all; close all
bounds = [ -55 55; %x
    -55 55; % y
    -pi pi]; % theta

dim = 3;
N_Nt = [2 3 4 5 6 7 8 9 10 11 12 13 14 15];

runs = 50;
rng('default');
sampleCountsActual = zeros(length(N_Nt)*runs,1);
for i = 1:length(N_Nt)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        samples = generatePertGridSamples(dim,N_Nt(i));
        sampleCountsActual(itr) = length(samples);

        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        samplesWrite = [samples(:,1:3)];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/PertGrid/3D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCountsActual(1))];
for i = 2:length(sampleCountsActual)
    str = [str,',',num2str(sampleCountsActual(i))];
end
disp([str '};']);
requested = N_Nt.^dim;
str = ['{' num2str(requested(1))];
for i = 2:length(requested)
    str = [str,',',num2str(requested(i))];
end
disp([str '};']);
disp('done');
%% To generate weighted rotated grid in OMPL for the SE(2) Maze
clear all; close all
bounds = [ -55 55; %x
    -55 55; % y
    -pi pi]; % theta

dim = 3;
N_Nt = [3 3 1;
    4 4 2;
    5 5 3;
    7 6 3;
    9 9 3;
    10 10 3;
    12 12 4;
    14 14 4;
    15 15 4;
    16 16 5;
    17 17 6;
    18 18 7;
    20 20 7;
    22 22 7];

runs = 50;
rng('default');
sampleCountsActual = zeros(length(N_Nt)*runs,1);
for i = 1:length(N_Nt)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        samples = generateWeightedRotatedGridSamples(N_Nt(i,:));
        sampleCountsActual(itr) = length(samples);

        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        samplesWrite = [samples(:,1:3)];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/SE2MazeRWGrid/3D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCountsActual(1))];
for i = 2:length(sampleCountsActual)
    str = [str,',',num2str(sampleCountsActual(i))];
end
disp([str '};']);
requested = N_Nt;
str = ['{' num2str(requested(1))];
for i = 2:length(requested)
    str = [str,',',num2str(prod(requested(i,:)))];
end
disp([str '};']);
disp('done');
%% To generate offset weighted grid in OMPL for the SE(2) Maze
clear all; close all
bounds = [ -55 55; %x
    -55 55; % y
    -pi pi]; % theta

dim = 3;
N_Nt = [3 3 1;
    4 4 2;
    5 5 3;
    7 6 3;
    9 9 3;
    10 10 3;
    12 12 4;
    14 14 4;
    15 15 4;
    16 16 5;
    17 17 6;
    18 18 7;
    20 20 7;
    22 22 7];

runs = 1;
rng('default');
sampleCountsActual = zeros(length(N_Nt)*runs,1);
for i = 1:length(N_Nt)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        samples = generateWeightedGridSamples(N_Nt(i,:));
        sampleCountsActual(itr) = length(samples);

        offset = [cos(2*pi*samples(:,3)) sin(2*pi*samples(:,3))];
        normInf = max(abs(offset),[],2);
        offset = offset./[normInf normInf];
        offset = 1/(3*N_Nt(i,1))*offset;
        samples(:,1:2) = samples(:,1:2) + offset;

        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        samplesWrite = [samples(:,1:3)];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/SE2MazeOffGrid_Hard/3D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCountsActual(1))];
for i = 2:length(sampleCountsActual)
    str = [str,',',num2str(sampleCountsActual(i))];
end
disp([str '};']);
requested = N_Nt;
str = ['{' num2str(requested(1))];
for i = 2:length(requested)
    str = [str,',',num2str(prod(requested(i,:)))];
end
disp([str '};']);
disp('done');
%% To generate offset weighted grid in OMPL for the SE(2) Maze
% Most current apr 24
clear all; close all
bounds = [ -55 55; %x
    -55 55; % y
    -pi pi]; % theta

dim = 3;
N_Nt = [3 1;
    4 2;
    5 2;
    6 2;
    7 2;
    8 2;
    9 2;
    8 3;
    9 3;
    10 3;
    12 3;
    14 3;
    15 3;
    16 3;
    17 3;
    18 3;
    20 3;
    22 3];

runs = 1;
rng('default');
sampleCountsActual = zeros(length(N_Nt)*runs,1);
for i = 1:length(N_Nt)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        samples = generateSE2Grid(dim,N_Nt(i,1),N_Nt(i,2));
        sampleCountsActual(itr) = length(samples);

        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        samplesWrite = [samples(:,1:3)];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/10piGrid/3D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCountsActual(1))];
for i = 2:length(sampleCountsActual)
    str = [str,',',num2str(sampleCountsActual(i))];
end
disp([str '};']);
requested = N_Nt(:,1).^2.*N_Nt(:,2).^2;
str = ['{' num2str(requested(1))];
for i = 2:length(requested)
    str = [str,',',num2str(requested(i,:))];
end
disp([str '};']);
disp('done');
%% Generate halton grid
clear all; close all
bounds = [ -55 55; %x
    -55 55; % y
    -pi pi]; % theta

dim = 3;
N_Nt = [3 1;
    4 2;
    5 2;
    6 2;
    7 2;
    8 2;
    9 2;
    8 3;
    9 3;
    10 3;
    12 3;
    14 3;
    15 3;
    16 3;
    17 3;
    18 3;
    20 3;
    22 3];

runs = 1;
rng('default');
sampleCountsActual = zeros(length(N_Nt)*runs,1);
for i = 1:length(N_Nt)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        samples = generateHaltonSamples(dim,N_Nt(i,1)^2*N_Nt(i,2)^2);
        sampleCountsActual(itr) = length(samples);

        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        samplesWrite = [samples(:,1:3)];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/Halton/3D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCountsActual(1))];
for i = 2:length(sampleCountsActual)
    str = [str,',',num2str(sampleCountsActual(i))];
end
disp([str '};']);
requested = N_Nt(:,1).^2.*N_Nt(:,2).^2;
str = ['{' num2str(requested(1))];
for i = 2:length(requested)
    str = [str,',',num2str(requested(i,:))];
end
disp([str '};']);
disp('done');
%% To generate rotated grid in OMPL for the SE(2) Maze
clear all; close all
bounds = [ -55 55; %x
    -55 55; % y
    -pi pi]; % theta

dim = 3;
N_Nt = [2 3 4 5 6 7 8 9 10 11 12 13 14 15];

runs = 1;
rng('default');
sampleCountsActual = zeros(length(N_Nt)*runs,1);
for i = 1:length(N_Nt)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        samples = generateRotatedGridSamples(dim,N_Nt(i));
        sampleCountsActual(itr) = length(samples);

        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        samplesWrite = [samples(:,1:3)];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/10piGrid/3D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCountsActual(1))];
for i = 2:length(sampleCountsActual)
    str = [str,',',num2str(sampleCountsActual(i))];
end
disp([str '};']);
requested = N_Nt.^dim;
str = ['{' num2str(requested(1))];
for i = 2:length(requested)
    str = [str,',',num2str(requested(i))];
end
disp([str '};']);
disp('done');