clearvars
clc
close all

% Customize Fonts
axisFontSize = 12;   % Font size for axis ticks
labelFontSize = 14;  % Font size for axis labels
legendFontSize = 12; % Font size for legend text
titleFontSize = 16;  % Font size for title


%% FIG 1: DUAL VALUE SENSITIVITY ANALYSIS PLOT %%

% Load data
DVs = table2array(readtable('DVs.csv'));
input = table2array(readtable('input.csv'));
rev = table2array(readtable('rev.csv'))/10e6; 
% Concatenate 2022 and 2023 together
DV5 = [DVs(:, 1); DVs(:, 2)];
DV1 = [DVs(:, 3); DVs(:, 4)];
DV13 = [DVs(:, 5); DVs(:, 6)];
DV2 = [DVs(:, 7); DVs(:, 8)];
DV3 = [DVs(:, 9); DVs(:, 10)];
DV4 = [DVs(:, 11); DVs(:, 12)];

LMP = [input(:, 1); input(:, 2)];
outflow = [input(:, 3); input(:, 4)];
inflow = [input(:, 5); input(:, 6)];
netflow = (outflow - inflow)/10e6;

% Create plot
% Create Figure with Specified Size (in inches)
figure('Units', 'inches', 'Position', [1, 1, 6, 4]);
hold on;
plot(DV5, 'DisplayName', '500', 'LineWidth', 1.5);
plot(DV1, 'DisplayName', '1000', 'LineWidth', 1.5);
plot(DV2, 'DisplayName', '2000', 'LineWidth', 1.5);
plot(DV3, 'DisplayName', '3000', 'LineWidth', 1.5);
plot(DV4, 'DisplayName', '4000', 'LineWidth', 1.5);
hold off;

% Add labels and legend
xlabel('Month','FontSize', labelFontSize,'FontName', 'Times New Roman');
ylabel('Price of Water ($/m^3)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');
xlim([0,24]);
xticks(0:3:24);
% title('2022-2023 Marginal Price of Water for Varying Feeder Capacity (FC)', 'FontSize', titleFontSize);
legend('FontSize', legendFontSize, 'Location', 'northeast','FontName', 'Times New Roman');
set(gca, 'FontSize', axisFontSize)
grid on;
saveas(gcf,'figures/waterprice_over_FC.png')

%% FIG 2: DUAL VALUE OVER LMP FOR FC = 1300 GW
figure('Units', 'inches', 'Position', [1, 1, 6, 4]);
yyaxis left;
plot(DV13, 'DisplayName', 'Price of Water', 'LineWidth', 1.5);
ylabel('Price of Water ($/m^3)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');

yyaxis right;
plot(LMP, 'DisplayName', 'Price of Electricity', 'LineWidth', 1.5);
ylabel('Price of Electricity ($/MWh)','FontSize', labelFontSize,'FontName', 'Times New Roman');

% Add labels and legend
xlabel('Month', 'FontSize', labelFontSize,'FontName', 'Times New Roman');
xlim([0,24]);
xticks(0:3:24);
% title('2022-2023 Marginal Price of Water over Marginal Price of Electricity', 'FontSize', titleFontSize);
legend('FontSize', legendFontSize-2.5, 'Location', 'northeast','FontName', 'Times New Roman');
set(gca, 'FontSize', axisFontSize)
grid on;
ax = gca; % Get the current axes
set(ax.YAxis, 'Color', 'k'); % Set the left y-axis ticks and labels to black
set(ax.YAxis(2), 'Color', 'k'); % Set the right y-axis ticks and labels to black

% Compute the Pearson correlation coefficient
r = corr(LMP, DV13); 
disp(r);

saveas(gcf,'figures/waterprice_over_LMP.png')

%% FIG 3: DUAL VALUE OVER NETFLOW FOR FC = 1300 GW
figure('Units', 'inches', 'Position', [1, 1, 6, 4]); 
yyaxis left;
plot(DV13, 'DisplayName', 'Price of Water', 'LineWidth', 1.5);
ylabel('Price of Water ($/m^3)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');

yyaxis right;
plot(netflow, 'DisplayName', 'Netflow', 'LineWidth', 1.5);
ylabel('Netflow (10^6 m^3)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');

% Add labels and legend
xlabel('Month', 'FontSize', labelFontSize,'FontName', 'Times New Roman');
xlim([0,24]);
xticks(0:3:24);
% title('2022-2023 Marginal Price of Water over System Netflow', 'FontSize', titleFontSize);
legend('FontSize', legendFontSize-2.5, 'Location', 'northeast','FontName', 'Times New Roman');
set(gca, 'FontSize', axisFontSize)
grid on;
saveas(gcf,'figures/waterprice_over_netflow.png')

%% FIG 4: Revenue over Feeder Capacity 
% Concatenate 2022 and 2023 together
rev5 = [rev(:, 1); rev(:, 2)];
rev1 = [rev(:, 3); rev(:, 4)];
rev13 = [rev(:, 5); rev(:, 6)];
rev2 = [rev(:, 7); rev(:, 8)];
rev3 = [rev(:, 9); rev(:, 10)];
rev4 = [rev(:, 11); rev(:, 12)];


% Create plot
figure('Units', 'inches', 'Position', [1, 1, 6, 4]); 
hold on;
plot(rev5, 'DisplayName', '500', 'LineWidth', 1.5);
plot(rev1, 'DisplayName', '1000', 'LineWidth', 1.5);
plot(rev2, 'DisplayName', '2000', 'LineWidth', 1.5);
plot(rev3, 'DisplayName', '3000', 'LineWidth', 1.5);
plot(rev4, 'DisplayName', '4000', 'LineWidth', 1.5);
hold off;

% Add labels and legend
xlabel('Month', 'FontSize', labelFontSize,'FontName', 'Times New Roman');
ylabel('Revenue (Million $)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');
xlim([0,24]);
xticks(0:3:24);
% title('2022-2023 Revenue for Varying Feeder Capacity (FC)', 'FontSize', titleFontSize);
legend('FontSize', legendFontSize, 'Location', 'northeast','FontName', 'Times New Roman');
set(gca, 'FontSize', axisFontSize)
grid on;
saveas(gcf,'figures/revenue_over_FC.png')