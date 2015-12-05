% Create Plots of simulation resoluts
%% Plot of convergence # points vs solution
clear all; close all;
% Data
optimal = sqrt(0.3^2+0.7^2) + sqrt(0.6^2+0.2^2);
randGrid = [...
    10	1.7293	1	0.0141
    20	1.6532	1	0.0157
    40	1.4974	1	0.054
    100	1.4419	1	0.1619
    200	1.4253	1	0.3506
    400	1.413	1	0.8117
    800	1.4082	1	2.0275
    1200    1.4052	1	3.8487
    2000    1.4016  1   9.9580];
randPoints = [...
    10	1.5909	0.275	0.0207
    20	1.5869	0.6	0.0254
    40	1.5722	0.875	0.0389
    100	1.4808	1	0.1248
    200	1.439	1	0.3036
    400	1.426	1	0.6699
    800	1.4113	1	1.7986
    1200	1.4069	1	3.6047
    2000    1.4038  1   9.4085];
gridPoints = [...
    10	1.7136	1	0.0579
    20	1.766	1	0.0132
    40	1.488	1	0.1113
    100	1.4944	1	0.1607
    200	1.4346	1	0.2097
    400	1.4137	1	0.8078
    800	1.4112	1	1.9324
    1200	1.4022	1	3.5421
    2000	1.4007	1	9.1247];
gridPoints16 = [...
    10	1.8216	0.85	0.0119
    20	1.8102	0.875	0.0154
    40	1.6126	1	0.0342
    100	1.5529	1	0.0973
    200	1.4742	1	0.2223
    400	1.4637	1	0.6447
    800	1.4447	1	1.9413
    1200	1.4424	1	3.8768];
gridPoints16noRand = [...
    10  1.7252  1   0.1903
    20  1.7839  1   0.0556
    40  1.5474  1   0.0966
    100	1.5662	1	0.0949
    200	1.5138	1	0.1871
    400	1.5363	1	0.6295
    800	1.5213	1	1.9647
    1200	1.5097	1	3.6303
    2000    1.5091  1   10.1208];

% Correct for actual number of samples
randGrid(:,1) = (floor(sqrt(randGrid(:,1)))+1).^2;
gridPoints(:,1) = (floor(sqrt(gridPoints(:,1)))+1).^2;

% Figure
figure; hold on;
plot([10 2000],[optimal optimal],'k--');
plot(randGrid(:,1),randGrid(:,2),'r-')
plot(randPoints(:,1),randPoints(:,2),'b-')
plot(gridPoints(:,1),gridPoints(:,2),'g-')
plot(gridPoints16(:,1),gridPoints16(:,2),'c-')
plot(gridPoints16noRand(:,1),gridPoints16noRand(:,2),'m-')

xlabel('Number of Samples');
ylabel('Cost');
legend('Optimal','Random Offset Grid','Random Sampling','Grid',...
    'Random Grid 16','Grid 16');
axis([0 1200 1.3 1.8]);

figure; hold on;
plot([10 2000],[optimal optimal],'k--');
plot(gridPoints16(:,1),gridPoints16(:,2),'c-')
plot(gridPoints16noRand(:,1),gridPoints16noRand(:,2),'m-')
legend('Optimal','Random Grid 16','Grid 16');
axis([0 1200 1.3 1.8]);