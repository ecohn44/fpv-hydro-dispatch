clearvars
clc
close all

% Customize Fonts
axisFontSize = 12;   % Font size for axis ticks
labelFontSize = 14;  % Font size for axis labels
legendFontSize = 12; % Font size for legend text
titleFontSize = 16;  % Font size for title

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


%% FIG 1: DUAL VALUE OVER LMP FOR FC = .5, 1, 2, 3 GW
figure('Units', 'inches', 'Position', [1, 1, 6, 4]);
yyaxis left;
%plot(DV13, 'DisplayName', 'Price of Water', 'LineWidth', 1.5);
hold on;
plot(DV5, 'b:', 'DisplayName', 'Price of Water, T = .5 GW', 'LineWidth', 3);
plot(DV1, 'b--', 'DisplayName', 'Price of Water, T = 1 GW', 'LineWidth', 2);
plot(DV2, 'b-', 'DisplayName', 'Price of Water, T = 2 GW', 'LineWidth', 1.5);
%plot(DV3, 'b-', 'DisplayName', 'Price of Water, FC = 3 GW', 'LineWidth', 1.5);
hold off;

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
r = corr(LMP, DV1); 
disp(r);

saveas(gcf,'figures/poster_fastwaterprice_over_LMP.png')


