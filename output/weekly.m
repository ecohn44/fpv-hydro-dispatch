clearvars
clc
close all

% Customize Fonts
axisFontSize = 12;   % Font size for axis ticks
labelFontSize = 14;  % Font size for axis labels
legendFontSize = 12; % Font size for legend text

FC = 1300 ;

% Read the CSV file
data = readtable('weekly_behavior.csv');

% Extract columns
p_s = data.x1;
p_h = data.x2;
u_t = data.x3;
p = p_h + p_s;

% Create a figure with four subplots
figure;
set(gcf, 'Position', [100, 100, 200, 800]);

% Plot p_h + p_s
subplot(2, 1, 1);
plot(p);
xlabel('Time (hour)');
xlim([0,168]);
ylabel('Total Generation (MW)');
yline(FC, 'r--', 'LineWidth', 1.5);

% Plot u_t
subplot(2, 1, 2);
plot(u_t);
xlabel('Time (hour)');
xlim([0,168]);
ylabel('Water Release Rate (m^3/s)', 'Interpreter', 'tex');
saveas(gcf,'figures/weekly.png')
