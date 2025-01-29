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
DVs = 100*table2array(readtable('fastDVs.csv')); %convert to cents
input = table2array(readtable('fastinput.csv'));
rev = table2array(readtable('fastrev.csv'))/10e6; %convert to M$
gen = table2array(readtable('fastgen.csv'));
s = gen(:,1);
h = gen(:,2);
combo = s + h;


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
ylabel('Price of Water (¢/m^3)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');
xlim([0,24]);
xticks(0:3:24);
% title('2022-2023 Marginal Price of Water for Varying Feeder Capacity (FC)', 'FontSize', titleFontSize);
legend('FontSize', legendFontSize, 'Location', 'northeast','FontName', 'Times New Roman');
set(gca, 'FontSize', axisFontSize)
grid on;
%saveas(gcf,'figures/fastwaterprice_over_FC.png')
savefig("figures/fastwaterprice_over_FC.fig")

%% FIG 2: DUAL VALUE OVER LMP FOR FC = 1300 GW
figure('Units', 'inches', 'Position', [1, 1, 6, 4]);
yyaxis left;
plot(DV13, 'DisplayName', 'Price of Water', 'LineWidth', 1.5);
ylabel('Price of Water (¢/m^3)','FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');

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

%saveas(gcf,'figures/fastwaterprice_over_LMP.png')
savefig("figures/fastwaterprice_over_LMP.fig")


%% FIG 3: Revenue over Feeder Capacity 
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
saveas(gcf,'figures/fastrevenue_over_FC.png')
savefig("figures/fastrevenue_over_FC.fig")



% Step 1: Sort the data
sorted_power = sort(combo);
sorted_hydro = sort(h);

% Step 2: Compute cumulative power
cumulative_power = normalize(cumsum(sorted_power), 'range', [0, 1]);
cumlative_hydro = normalize(cumsum(sorted_hydro), 'range', [0, 1]);

% Step 3: Normalize x-axis (percentage of time)
time_percentage = linspace(0, 100, length(sorted_power));

% Step 4: Plot the CDF
figure;
hold on 
plot(time_percentage, cumulative_power, 'LineWidth', 2, 'DisplayName', 'Total Power');
plot(time_percentage, cumlative_hydro, 'LineWidth', 2, 'DisplayName', 'Hydropower');
xlabel('Percentage of Time (%)');
ylabel('Cumulative Generation (normalized)');
title('CDF of Power Generation');
legend('FontSize', legendFontSize, 'Location', 'northwest','FontName', 'Times New Roman');


% Step 5: empirical CDF
figure;
hold on
ecdf(combo);
hold on
ecdf(h);
legend('Combo', 'Hydropower', 'Location', 'best');
xlabel('Power Generation (MW)');
ylabel('Empirical CDF');
title('Empirical CDF of Time Series Data');

%{ Interpretation (Dr Lall)
% 50 percent of time hydro is 600MW and combo is 1300 MW
% I am used to plotting this with Power on the y axis and time on the x axis
% In this case it looks like combo is 1300 with cdf =0.3 , i.e 30 of time you produces less than 1300
% you should be plotting solar and combo since the point is to show how the reliability of solr is improved
%}




%{
%% DECOMM: DUAL VALUE OVER NETFLOW FOR FC = 1300 GW
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
%}