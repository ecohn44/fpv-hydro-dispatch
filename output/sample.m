% Sample Data
x = linspace(0, 10, 100);
y1 = sin(x);
y2 = cos(x);

% Create Figure with Specified Size (in inches)
figure('Units', 'inches', 'Position', [1, 1, 6, 4]);
% Plot Data
plot(x, y1, '-b', 'LineWidth', 1.5);
hold on;
plot(x, y2, '--r', 'LineWidth', 1.5);
hold off;

% Customize Fonts
axisFontSize = 12;   % Font size for axis ticks
labelFontSize = 14;  % Font size for axis labels
legendFontSize = 12; % Font size for legend text
titleFontSize = 16;  % Font size for title

% Add Labels and Title
xlabel('X Axis', 'FontSize', labelFontSize);
ylabel('Y Axis', 'FontSize', labelFontSize);
title('Sample Plot with Controlled Font Sizes', 'FontSize', titleFontSize);

% Add Legend
legend({'sin(x)', 'cos(x)'}, 'FontSize', legendFontSize, 'Location', 'best');
% Set Axis Font Size
set(gca, 'FontSize', axisFontSize);
% Display Grid
grid on;
% Show the plot
