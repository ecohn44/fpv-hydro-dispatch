clearvars
clc
close all

% Customize Fonts
axisFontSize = 12;   % Font size for axis ticks
labelFontSize = 14;  % Font size for axis labels
legendFontSize = 12; % Font size for legend text

FC = 1300 ;

% Read the CSV file
data = readtable('relaxed_weekly_behavior.csv');
b_data = readtable('baseline_weekly_behavior.csv');

% Extract columns
p_s = data.x1;
p_h = data.x2;
u_t = data.x3;
p = p_h + p_s;

% Extract columns
bp_s = b_data.x1;
bp_h = b_data.x2;
bu_t = b_data.x3;
bp = bp_h + bp_s;

figure('Units', 'inches', 'Position', [1, 1, 6, 8]); % Double height for two subplots stacked vertically

%{
% Plot p_h + p_s
subplot(2, 1, 1);
hold on;
plot(p, '-o', 'DisplayName', 'Relaxed Generation', 'LineWidth', 1.5,'MarkerFaceColor', 'b','MarkerSize', 3);
plot(bp,'-.', 'DisplayName', 'Baseline Generation', 'LineWidth', 1.5);
hold off;
xlabel('Time (hour)','FontName', 'Times New Roman','FontSize', labelFontSize);
xlim([0,168]);
xticks(0:24:168);
ylabel('Total Generation (MW)','FontName', 'Times New Roman','FontSize', labelFontSize);
yline(FC, 'k--', 'DisplayName', 'FC','LineWidth', 2.0);
set(gca, 'FontSize', axisFontSize)
%set(gca, 'Units', 'inches', 'Position', [0.75, 4.5, 5.5, 3.5]); % Position to fit within figure
%set(gcf, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman');
%legend('FontSize', legendFontSize, 'Location', 'southwest','FontName', 'Times New Roman');
%}

% Plot p_s
subplot(3, 1, 1);
hold on;
plot(bp_s,'-o', 'DisplayName', 'M1', 'LineWidth', 1.5,'MarkerFaceColor', 'b','MarkerSize', 3);
plot(p_s, '-.', 'DisplayName', 'M2', 'LineWidth', 1.5);
hold off;
xlabel('Time (h)','FontName', 'Times New Roman','FontSize', labelFontSize);
xlim([0,168]);
ylim([0,1500]);
xticks(0:24:168);
ylabel('FPV Generation (MW)','FontName', 'Times New Roman','FontSize', labelFontSize);
legend('FontSize', legendFontSize, 'Location', 'northeast','FontName', 'Times New Roman','Box','off');
set(gca, 'FontSize', axisFontSize)
%set(gca, 'Units', 'inches', 'Position', [0.75, 4.5, 5.5, 3.5]); % Position to fit within figure
%set(gcf, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman');
%legend('FontSize', legendFontSize, 'Location', 'southwest','FontName', 'Times New Roman');

% plot p_h
subplot(3, 1, 2);
hold on;
plot(bp_h, 'DisplayName', 'M1', 'LineWidth', 1.5,'MarkerFaceColor', 'b','MarkerSize', 3);
plot(p_h, '-.', 'DisplayName', 'M2', 'LineWidth', 1.5);
hold off;
xlabel('Time (h)','FontName', 'Times New Roman','FontSize', labelFontSize);
xlim([0,168]);
%ylim([0,1500]);
xticks(0:24:168);
ylabel('Hydropower Generation (MW)','FontName', 'Times New Roman','FontSize', labelFontSize);
legend('FontSize', legendFontSize, 'Location', 'northeast','FontName', 'Times New Roman','Box','off');
set(gca, 'FontSize', axisFontSize)

% Plot u_t
subplot(3, 1, 3);
hold on
plot(bu_t, 'DisplayName', 'M1', 'LineWidth', 1.5,'MarkerFaceColor', 'b','MarkerSize', 3);
plot(u_t, '-.', 'DisplayName', 'M2', 'LineWidth', 1.5);
hold off
xlabel('Time (h)','FontSize', labelFontSize,'FontName', 'Times New Roman');
xlim([0,168]);
ylim([0, 600])
xticks(0:24:168);
ylabel('Water Release Rate (m^3/s)', 'Interpreter', 'tex','FontName', 'Times New Roman','FontSize', labelFontSize);
legend('FontSize', legendFontSize, 'Location', 'northeast','FontName', 'Times New Roman','Box','off');
set(gca, 'FontSize', axisFontSize)
%set(gca, 'Units', 'inches', 'Position', [0.75, 0.5, 5.5, 3.5]); % Position to fit within figure
%set(gcf, 'DefaultAxesFontName', 'Times New Roman', 'DefaultTextFontName', 'Times New Roman');
%legend('FontSize', legendFontSize, 'Location', 'southwest','FontName', 'Times New Roman');
saveas(gcf,'figures/weekly.png')
savefig("figures/weekly.fig")
