clearvars
clc
close all

%% DUAL VALUE SENSITIVITY ANALYSIS PLOT %%

% Load data
DVs = table2array(readtable('DVs.csv'));
input = table2array(readtable('input.csv'));

% Concatenate 2022 and 2023 together
DV1 = [DVs(:, 1); DVs(:, 2)];
DV2 = [DVs(:, 3); DVs(:, 4)];
DV3 = [DVs(:, 5); DVs(:, 6)];
DV4 = [DVs(:, 7); DVs(:, 8)];
LMP = [input(:, 1); input(:, 2)];
outflow = [input(:, 3); input(:, 4)];
inflow = [input(:, 5); input(:, 6)];
netflow = (outflow - inflow)/10e6;

% Create plot
figure;
hold on;
plot(DV1, 'DisplayName', 'PF = 1000', 'LineWidth', 1.5);
plot(DV2, 'DisplayName', 'PF = 2000', 'LineWidth', 1.5);
plot(DV3, 'DisplayName', 'PF = 3000', 'LineWidth', 1.5);
plot(DV4, 'DisplayName', 'PF = 4000', 'LineWidth', 1.5);
hold off;

% Add labels and legend
xlabel('Month');
ylabel('Price of Water $/m3');
xlim([0,24]);
title('2022-2023 Marginal Price of Water for Varying Feeder Capacity (PF)', 'FontSize', 16);
legend show;
grid on;

%% DUAL VALUE OVER LMP FOR PF = 3000 GW
figure;
yyaxis left;
plot(DV3, 'DisplayName', 'Price of Water ($/m3)', 'LineWidth', 1.5);
ylabel('Price of Water $/m3');

yyaxis right;
plot(LMP, 'DisplayName', 'Price of Electricity ($/MW)', 'LineWidth', 1.5);
ylabel('Price of Electricity $/MW');

% Add labels and legend
xlabel('Month');
xlim([0,24]);
title('2022-2023 Marginal Price of Water over Marginal Price of Electricity', 'FontSize', 16);
legend show;
grid on;

%% DUAL VALUE OVER NETFLOW FOR PF = 3000 GW
figure;
yyaxis left;
plot(DV3, 'DisplayName', 'Price of Water ($/m3)', 'LineWidth', 1.5);
ylabel('Price of Water $/m3');

yyaxis right;
plot(netflow, 'DisplayName', 'Netflow (Mm3)', 'LineWidth', 1.5);
ylabel('Netflow (Million m3)');

% Add labels and legend
xlabel('Month');
xlim([0,24]);
title('2022-2023 Marginal Price of Water over System Netflow', 'FontSize', 16);
legend show;
grid on;