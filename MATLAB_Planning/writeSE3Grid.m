%% Full grid (weighted, rotated)
clear all; close all;
% cubicles full
bounds = [319.62  -508.88; % x
	531.87 -230.13; % y 
	101.00 -123.75]; % z 
% cubicles short
% bounds = [319.62  -508.88; % x
% 	531.87 -230.13; % y
% 	101.00 0]; % z

dim = 6;
sides = [   ...
    2 2 2 1 1 1 % 8
    2 2 2 2 2 1 % 32
    2 2 2 2 2 2 % 64
    3 3 2 2 2 2 % 144
    3 3 3 2 2 2 % 216
    3 3 3 3 2 2 % 324
    3 3 3 3 3 2 % 486
    3 3 3 3 3 3 % 729
    4 3 3 3 3 3 % 972
    4 4 3 3 3 3 % 1296
    4 4 4 3 3 3 % 1728
    4 4 4 4 3 3 % 2304
    4 4 4 4 4 3 % 3072
    4 4 4 4 4 4 % 4096
    5 4 4 4 4 4 % 5120
    5 5 4 4 4 4]; % 6400
%     5 5 5 4 4 4 % 8000
%     5 5 5 5 4 4 % 10000
%     5 5 5 5 5 4 % 12500
%     5 5 5 5 5 5]; % 15625
runs = 1;
rng('default');
sampleCounts = zeros(length(sides)*runs,1);
for i = 1:length(sides)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        side = sides(i,:);
        samples = generateWeightedRotatedGridSamples(side);
        sampleCounts(itr) = length(samples);
        phi = 2*pi*samples(:,4);
        theta = acos(2*samples(:,5)-1);
        psi = pi*(2*samples(:,6)-1);
        
        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        % angle-axis
        ax = sin(phi/2).*sin(theta/2).*cos(psi/2) + cos(phi/2).*cos(theta/2).*sin(psi/2);
        ay = sin(phi/2).*cos(theta/2).*cos(psi/2) + cos(phi/2).*sin(theta/2).*sin(psi/2);
        az = cos(phi/2).*sin(theta/2).*cos(psi/2) - sin(phi/2).*cos(theta/2).*sin(psi/2);
        normfactor = sqrt(ax.^2 + ay.^2 + az.^2);
        ax = ax./normfactor;
        ay = ay./normfactor;
        az = az./normfactor;
        angles = 2*acos(cos(phi/2).*cos(theta/2).*cos(psi/2) - sin(phi/2).*sin(theta/2).*sin(psi/2));
        
        % quaternions
        q = zeros(length(phi),4);
        for j = 1:length(phi)
            [q(j,1), q(j,2), q(j,3), q(j,4)] = ...
                euler2quat(phi(j),theta(j),psi(j),'xyz');
        end
        
        samplesWrite = [samples(:,1:3) ax ay az angles];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/10piGrid/7D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCounts(1))];
for i = 2:length(sampleCounts)
    str = [str,',',num2str(sampleCounts(i))];
end
disp([str '};']);
disp('done');
%% Perturbed Grid (weighted)
clear all; close all;
% cubicles full
bounds = [319.62  -508.88; % x
    531.87 -230.13; % y 
    101.00 -123.75]; % z 
% cubicles short
% bounds = [319.62  -508.88; % x
%     531.87 -230.13; % y
%     101.00 0]; % z

dim = 6;
sides = [   ...
    2 2 2 1 1 1 % 8
    2 2 2 2 2 1 % 32
    2 2 2 2 2 2 % 64
    3 3 2 2 2 2 % 144
    3 3 3 2 2 2 % 216
    3 3 3 3 2 2 % 324
    3 3 3 3 3 2 % 486
    3 3 3 3 3 3 % 729
    4 3 3 3 3 3 % 972
    4 4 3 3 3 3 % 1296
    4 4 4 3 3 3 % 1728
    4 4 4 4 3 3 % 2304
    4 4 4 4 4 3 % 3072
    4 4 4 4 4 4 % 4096
    5 4 4 4 4 4 % 5120
    5 5 4 4 4 4]; % 6400
%     5 5 5 4 4 4 % 8000
%     5 5 5 5 4 4 % 10000
%     5 5 5 5 5 4 % 12500
%     5 5 5 5 5 5]; % 15625
runs = 50;
rng('default');
sampleCounts = zeros(length(sides)*runs,1);
for i = 1:length(sides)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        side = sides(i,:);
        samples = generateWeightedPertGridSamples(side);
        sampleCounts(itr) = length(samples);
        phi = 2*pi*samples(:,4);
        theta = asin(samples(:,5));
        psi = pi*(2*samples(:,6)-1);
        
        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        % angle-axis
        ax = sin(phi/2).*sin(theta/2).*cos(psi/2) + cos(phi/2).*cos(theta/2).*sin(psi/2);
        ay = sin(phi/2).*cos(theta/2).*cos(psi/2) + cos(phi/2).*sin(theta/2).*sin(psi/2);
        az = cos(phi/2).*sin(theta/2).*cos(psi/2) - sin(phi/2).*cos(theta/2).*sin(psi/2);
        normfactor = sqrt(ax.^2 + ay.^2 + az.^2);
        ax = ax./normfactor;
        ay = ay./normfactor;
        az = az./normfactor;
        angles = 2*acos(cos(phi/2).*cos(theta/2).*cos(psi/2) - sin(phi/2).*sin(theta/2).*sin(psi/2));
        
        % quaternions
        q = zeros(length(phi),4);
        for j = 1:length(phi)
            [q(j,1), q(j,2), q(j,3), q(j,4)] = ...
                euler2quat(phi(j),theta(j),psi(j),'xyz');
        end
        
        samplesWrite = [samples(:,1:3) ax ay az angles];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/PertGrid/7D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCounts(1))];
for i = 2:length(sampleCounts)
    str = [str,',',num2str(sampleCounts(i))];
end
disp([str '};']);
disp('done');
%% Full grid (weighted, rotated, spread) *most current apr 23 @ 9pm
clear all; close all;
% cubicles full
bounds = [319.62  -508.88; % x
	531.87 -230.13; % y 
	101.00 -123.75]; % z 
% cubicles short
% bounds = [319.62  -508.88; % x
% 	531.87 -230.13; % y
% 	101.00 0]; % z

dim = 6;
sides = [   ...
    2 2 2 1 1 1 % 8
    2 2 2 2 2 1 % 32
    2 2 2 2 2 2 % 64
    3 3 2 2 2 2 % 144
    3 3 3 2 2 2 % 216
    3 3 3 3 2 2 % 324
    3 3 3 3 3 2 % 486
    3 3 3 3 3 3 % 729
    4 3 3 3 3 3 % 972
    4 4 3 3 3 3 % 1296
    4 4 4 3 3 3 % 1728
    4 4 4 4 3 3 % 2304
    4 4 4 4 4 3 % 3072
    4 4 4 4 4 4 % 4096
    5 4 4 4 4 4 % 5120
    5 5 4 4 4 4]; % 6400
%     5 5 5 4 4 4 % 8000
%     5 5 5 5 4 4 % 10000
%     5 5 5 5 5 4 % 12500
%     5 5 5 5 5 5]; % 15625
runs = 1;
rng('default');
sampleCounts = zeros(length(sides)*runs,1);
for i = 1:length(sides)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        side = sides(i,:);
        samples = generateSE3Grid(side);
        
        sampleCounts(itr) = size(samples,1);
        phi = 2*pi*samples(:,4);
        theta = asin(samples(:,5));
        psi = pi*(2*samples(:,6)-1);
        
        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        % angle-axis
        ax = sin(phi/2).*sin(theta/2).*cos(psi/2) + cos(phi/2).*cos(theta/2).*sin(psi/2);
        ay = sin(phi/2).*cos(theta/2).*cos(psi/2) + cos(phi/2).*sin(theta/2).*sin(psi/2);
        az = cos(phi/2).*sin(theta/2).*cos(psi/2) - sin(phi/2).*cos(theta/2).*sin(psi/2);
        normfactor = sqrt(ax.^2 + ay.^2 + az.^2);
        ax = ax./normfactor;
        ay = ay./normfactor;
        az = az./normfactor;
        angles = 2*acos(cos(phi/2).*cos(theta/2).*cos(psi/2) - sin(phi/2).*sin(theta/2).*sin(psi/2));
        
        % quaternions
        q = zeros(length(phi),4);
        for j = 1:length(phi)
            [q(j,1), q(j,2), q(j,3), q(j,4)] = ...
                euler2quat(phi(j),theta(j),psi(j),'xyz');
        end
        
        samplesWrite = [samples(:,1:3) ax ay az angles];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/10piGrid/7D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCounts(1))];
for i = 2:length(sampleCounts)
    str = [str,',',num2str(sampleCounts(i))];
end
disp([str '};']);
disp('done');
%% Halton Sampling
clear all; close all;
% cubicles full
bounds = [319.62  -508.88; % x
	531.87 -230.13; % y 
	101.00 -123.75]; % z 
% cubicles short
% bounds = [319.62  -508.88; % x
% 	531.87 -230.13; % y
% 	101.00 0]; % z

dim = 6;
sides = [   ...
    2 2 2 1 1 1 % 8
    2 2 2 2 2 1 % 32
    2 2 2 2 2 2 % 64
    3 3 2 2 2 2 % 144
    3 3 3 2 2 2 % 216
    3 3 3 3 2 2 % 324
    3 3 3 3 3 2 % 486
    3 3 3 3 3 3 % 729
    4 3 3 3 3 3 % 972
    4 4 3 3 3 3 % 1296
    4 4 4 3 3 3 % 1728
    4 4 4 4 3 3 % 2304
    4 4 4 4 4 3 % 3072
    4 4 4 4 4 4 % 4096
    5 4 4 4 4 4 % 5120
    5 5 4 4 4 4]; % 6400
%     5 5 5 4 4 4 % 8000
%     5 5 5 5 4 4 % 10000
%     5 5 5 5 5 4 % 12500
%     5 5 5 5 5 5]; % 15625
runs = 1;
rng('default');
sampleCounts = zeros(length(sides)*runs,1);
for i = 1:length(sides)
    for k = 1:runs
        rng(k);
        itr = (i-1)*runs+k
        side = sides(i,:);
        samples = generateHaltonSamples(dim,prod(side));
        
        sampleCounts(itr) = size(samples,1);
        phi = 2*pi*samples(:,4);
        theta = asin(samples(:,5));
        psi = pi*(2*samples(:,6)-1);
        
        for j = 1:size(bounds,1)
            high = bounds(j,1);
            low = bounds(j,2);
            delta = high-low;
            samples(:,j) = delta*samples(:,j) + low;
        end
        
        % angle-axis
        ax = sin(phi/2).*sin(theta/2).*cos(psi/2) + cos(phi/2).*cos(theta/2).*sin(psi/2);
        ay = sin(phi/2).*cos(theta/2).*cos(psi/2) + cos(phi/2).*sin(theta/2).*sin(psi/2);
        az = cos(phi/2).*sin(theta/2).*cos(psi/2) - sin(phi/2).*cos(theta/2).*sin(psi/2);
        normfactor = sqrt(ax.^2 + ay.^2 + az.^2);
        ax = ax./normfactor;
        ay = ay./normfactor;
        az = az./normfactor;
        angles = 2*acos(cos(phi/2).*cos(theta/2).*cos(psi/2) - sin(phi/2).*sin(theta/2).*sin(psi/2));
        
        % quaternions
        q = zeros(length(phi),4);
        for j = 1:length(phi)
            [q(j,1), q(j,2), q(j,3), q(j,4)] = ...
                euler2quat(phi(j),theta(j),psi(j),'xyz');
        end
        
        samplesWrite = [samples(:,1:3) ax ay az angles];
        filename = ['/Users/brianichter/Desktop/ASL_Shared/OMPL_Results/ISRR15/Halton/7D/48/Samples_',num2str(itr-1),'.txt'];
        dlmwrite(filename,dim);
        dlmwrite(filename,samplesWrite,'-append','delimiter',' ');
    end
end

str = ['{' num2str(sampleCounts(1))];
for i = 2:length(sampleCounts)
    str = [str,',',num2str(sampleCounts(i))];
end
disp([str '};']);
disp('done');