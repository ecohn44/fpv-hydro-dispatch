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
figure('Units', 'inches', 'Position', [1, 1, 6, 8]); % Double height for two subplots stacked vertically

% Plot p_h + p_s
subplot(2, 1, 1);
plot(p,'DisplayName', 'Generation', 'LineWidth', 1.5);
xlabel('Time (hour)','FontName', 'Times New Roman','FontSize', labelFontSize);
xlim([0,168]);
xticks(0:24:168);
ylabel('Total Generation (MW)','FontName', 'Times New Roman','FontSize', labelFontSize);
yline(FC, 'r--', 'DisplayName', 'FC','LineWidth', 1.5);
grid on;
set(gca, 'FontSize', axisFontSize)
%set(gca, 'Units', 'inches', 'Position', [0.75, 4.5, 5.5, 3.5]); % Position to fit within figure
%set(gcf, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman');
%legend('FontSize', legendFontSize, 'Location', 'southwest','FontName', 'Times New Roman');

% Plot u_t
subplot(2, 1, 2);
plot(u_t,'DisplayName', 'Water Release', 'LineWidth', 1.5);
xlabel('Time (hour)','FontSize', labelFontSize,'FontName', 'Times New Roman');
xlim([0,168]);
xticks(0:24:168);
ylabel('Water Release Rate (m^3/s)', 'Interpreter', 'tex','FontName', 'Times New Roman','FontSize', labelFontSize);
grid on;
set(gca, 'FontSize', axisFontSize)
%set(gca, 'Units', 'inches', 'Position', [0.75, 0.5, 5.5, 3.5]); % Position to fit within figure
%set(gcf, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman');
%legend('FontSize', legendFontSize, 'Location', 'southwest','FontName', 'Times New Roman');
saveas(gcf,'figures/weekly.png')
