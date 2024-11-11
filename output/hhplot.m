%% Lake Mead Data Import

% Customize Fonts
axisFontSize = 12;   % Font size for axis ticks
labelFontSize = 14;  % Font size for axis labels
legendFontSize = 12; % Font size for legend text
titleFontSize = 16;  % Font size for title

% Total Capacity of Lake Mead (Ac-Ft)
mead_vol = [30167000; 28946000; 28667000; 26483000; 10230000; 4553000; 2547000];
meadlogvol = log(mead_vol);
mead_vol_m = mead_vol*1233.48; % convert Ac-Ft to cubic meters

% Elevation of Lake Mead (Ft)
mead_h = [1229.0; 1221.4; 1219.6; 1205.4; 1050.0; 950.0; 895.0];
meadlogh = log(mead_h);
mead_h_m = mead_h*0.3048; % convert ft to m

a = 14.9837;
b = 0.1321;
% Generate fitted line data
x_fit = linspace(min(mead_vol_m), max(mead_vol_m), 100); % Generate x values for smooth line
y_fit = a * x_fit.^b; % Compute corresponding y values



% Create plot
figure('Units', 'inches', 'Position', [1, 1, 6, 4]); 
hold on;
scatter(mead_vol_m/10e9, mead_h_m, 100, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Input Data');
plot(x_fit/10e9, y_fit, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Best Fit'); % Red line for best fit
hold off;

% Add labels and legend
xlabel('Volume (10^9 m^3)', 'FontSize', labelFontSize, 'Interpreter', 'tex','FontName', 'Times New Roman');
ylabel('Height (m)','FontSize', labelFontSize,'FontName', 'Times New Roman');
legend('FontSize', legendFontSize, 'Location', 'northwest','FontName', 'Times New Roman');
set(gca, 'FontSize', axisFontSize)
grid on;
saveas(gcf,'figures/hh_fit.png')
savefig("figures/hh_fit.fig")