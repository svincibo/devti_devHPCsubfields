% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields/';

% Select WM measure.
wm = 'md';
subregion = {'b_ca1', 'b_ca23', 'b_sub'};
degp = 1;

% Bring in Meg's data.
data_meg = readtable(fullfile(rootDir, 'supportFiles/dti_subfields_long.csv'));
sub_include = unique(data_meg.subnr);

% Bring in my data.
load(fullfile(rootDir, 'supportFiles/devti_data_md_mrtrix3act.mat')); m = m*1000; 
% For fits from FSLDTIFIT, use: devti_data_md.mat
% For fits from mrtrix3 act, use: devti_data_md_mrtrix3act.mat %NOTE: need to scale m by 1e3 on line 21 to match.

% Convert data to table for easier model specification.
data = array2table(cat(2, transpose(sub), transpose(age), transpose(sex), transpose(iq), m*1000), 'VariableNames', {'subID', 'age', 'sex',  'iq', roi{1, :}});

% Scale Meg's data up because brainlife.io outputs it in 10e-3.
data_meg.md = data_meg.md*1000;

keep_mydata = find(ismember(data.subID, data_meg.subnr));
keep_megdata_ca1 = find(ismember(data_meg.subnr, data.subID) & strcmp(data_meg.subfield, 'ca1'));
keep_megdata_ca23 = find(ismember(data_meg.subnr, data.subID) & strcmp(data_meg.subfield, 'ca23'));
keep_megdata_sub = find(ismember(data_meg.subnr, data.subID) & strcmp(data_meg.subfield, 'sub'));

%% Visualize: Age x Meg's Data

figure(1)
hold on;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
alphablend = .8;
ylim_lo = 0.5;
ylim_hi = 1.2;
ytickvalues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
xlim_lo = 0;
xlim_hi = 30;
xtickvalues = [xlim_lo xlim_hi];

% Correlations between age and average WM measure in ROI.
% CA1
x = data_meg.age(find(strcmp(data_meg.subfield, 'ca1')));
y = data_meg.md(find(strcmp(data_meg.subfield, 'ca1')));
c1_color = [0 0.4470 0.7410];
scatter(x, y, 'filled', 'MarkerEdgeColor', c1_color, 'MarkerFaceColor', c1_color, 'SizeData', markersize)

% c1 = corr(x(~isnan(y)), y(~isnan(y)));
c1 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
x1 = linspace(0,30);
y1 = 1./(1+x1);

f1 = polyval(c1,x1);
plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c1_color)
clear x
hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c1_color(1:3));
set(hp3, 'facecolor', c1_color(1:3), 'edgecolor', 'none', 'facealpha', .2);

disp(['CA1 = ' num2str(c1)])
clear x y f1 hp3

% CA23dg
x = data_meg.age(find(strcmp(data_meg.subfield, 'ca23')));
y = data_meg.md(find(strcmp(data_meg.subfield, 'ca23')));
c2_color = [0.4940 0.1840 0.5560];
scatter(x(~isnan(y)), y(~isnan(y)), 'filled', 'MarkerEdgeColor', c2_color, 'MarkerFaceColor', c2_color, 'SizeData', markersize)

% c2 = corr(x(~isnan(y)), y(~isnan(y)));
c2 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);

f1 = polyval(c2,x1);
plot(x1,f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c2_color)
clear x
hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c2_color(1:3));
set(hp3, 'facecolor', c2_color(1:3), 'edgecolor', 'none', 'facealpha', .2);

disp(['CA23 = ' num2str(c2)])
clear x y f1 hp3

% SUB
x = data_meg.age(find(strcmp(data_meg.subfield, 'sub')));
y = data_meg.md(find(strcmp(data_meg.subfield, 'sub')));
c3_color = [0.8500 0.3250 0.0980];
scatter(x(~isnan(y)), y(~isnan(y)), 'filled', 'MarkerEdgeColor', c3_color, 'MarkerFaceColor', c3_color, 'SizeData', markersize)

% c3 = corr(x(~isnan(y)), y(~isnan(y)));
c3 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);

f1 = polyval(c3,x1);
plot(x1,f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c3_color)
clear x
hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c3_color(1:3));
set(hp3, 'facecolor', c3_color(1:3), 'edgecolor', 'none', 'facealpha', .2);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
% xax.Limits = [5 30];
% xax.TickValues = [10 15 20 25 30];
% xax.TickLabels = {'10', '15', '20', '25', '30'};

xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
% yax.TickValues = [ylim_lo .75 1];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.1f'), '', num2str(ylim_hi, '%1.1f')};
% yax.TickLabels = {'0.5', '0.75', '1'};

yax.FontName = fontname;
yax.FontSize = fontsize;

legend([{'CA1'}, {''}, {''}, {'CA23'}, {''}, {''}, {'SUB'}, {''}, {''}], 'Location', 'southwest')
legend box off

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

% a.YLabel.String = {'Meg''s Data'; '(MD x 1e-3)'};
a.YLabel.String = {'Mean Diffusivity (MD), (MD x 1e-3)'};

a.XLabel.String = {'Age'; '(years)'};
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_equality_md_meg_age'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_equality_md_meg_age'), '-depsc')

hold off;

%% Visualize:Age x My Data FSLDTIFIT

figure(2)
hold on;

% Correlations between age and average WM measure in ROI.
% CA1
x = transpose(age(keep_mydata));
y = m(keep_mydata, strcmp(roi(1, :), 'b_ca1'));
c1_color = [0 0.4470 0.7410];
scatter(x, y, 'filled', 'MarkerEdgeColor', c1_color, 'MarkerFaceColor', c1_color, 'SizeData', markersize)

% c1 = corr(x(~isnan(y)), y(~isnan(y)));
c1 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
x1 = linspace(0,30);
y1 = 1./(1+x1);

f1 = polyval(c1,x1);
plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c1_color)
clear x
hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c1_color(1:3));
set(hp3, 'facecolor', c1_color(1:3), 'edgecolor', 'none', 'facealpha', .2);

disp(['CA1 = ' num2str(c1)])
clear x y f1 hp3

% CA23dg
x = transpose(age(keep_mydata));
y = m(keep_mydata, strcmp(roi(1, :), 'b_ca23dg'));
c2_color = [0.4940 0.1840 0.5560];
scatter(x(~isnan(y)), y(~isnan(y)), 'filled', 'MarkerEdgeColor', c2_color, 'MarkerFaceColor', c2_color, 'SizeData', markersize)

% c2 = corr(x(~isnan(y)), y(~isnan(y)));
c2 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);

f1 = polyval(c2,x1);
plot(x1,f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c2_color)
clear x
hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c2_color(1:3));
set(hp3, 'facecolor', c2_color(1:3), 'edgecolor', 'none', 'facealpha', .2);

disp(['CA23 = ' num2str(c2)])
clear x y f1 hp3

% SUB
x = transpose(age(keep_mydata));
y = m(keep_mydata, strcmp(roi(1, :), 'b_sub'));
c3_color = [0.8500 0.3250 0.0980];
scatter(x(~isnan(y)), y(~isnan(y)), 'filled', 'MarkerEdgeColor', c3_color, 'MarkerFaceColor', c3_color, 'SizeData', markersize)

% c3 = corr(x(~isnan(y)), y(~isnan(y)));
c3 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);

f1 = polyval(c3,x1);
plot(x1,f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c3_color)
clear x
hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c3_color(1:3));
set(hp3, 'facecolor', c3_color(1:3), 'edgecolor', 'none', 'facealpha', .2);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.1f'), '', num2str(ylim_hi, '%1.1f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

legend([{'CA1'}, {''}, {''}, {'CA23'}, {''}, {''}, {'SUB'}, {''}, {''}], 'Location', 'southwest')
legend box off

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

% a.YLabel.String = {'Sophia''s Data'; '(MD x 1e-3)'};
a.YLabel.String = {'Mean Diffusivity (MD), (MD x 1e-3)'};

a.XLabel.String = {'Age'; '(years)'};
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_equality_md_bl_age_mrtrix3act'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_equality_md_bl_age_mrtrix3act'), '-depsc')

hold off;
