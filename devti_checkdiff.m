% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields/';

% Bring in Meg's data.
data_meg = readtable(fullfile(rootDir, 'supportFiles/dti_subfields_long.csv'));
sub_include = unique(data_meg.subnr);

% Bring in my data.
load(fullfile(rootDir, 'supportFiles/devti_data_md_mrtrix3act.mat'))
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

%% Visualize.

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
xylim_lo = 0.5;
xylim_hi = 1.4;
xtickvalues = [xylim_lo (xylim_lo+xylim_hi)/2 xylim_hi];

plot(linspace(xylim_lo, xylim_hi, 2), linspace(xylim_lo, xylim_hi, 2), ':k')
hold on;

% Get md diff for ca1.
x = data.b_ca1(keep_mydata);
y = data_meg.md(keep_megdata_ca1);
c1_color = [0 0.4470 0.7410];
plot(x, y, 'o', 'MarkerEdgeColor', c1_color, 'MarkerFaceColor', c1_color)

% Get md difffor ca23.
x = data.b_ca23(keep_mydata);
y = data_meg.md(keep_megdata_ca23);
c2_color = [0.4940 0.1840 0.5560];
plot(x, y, 'o', 'MarkerEdgeColor', c2_color, 'MarkerFaceColor', c2_color)

% Get md difffor sub.
x = data.b_sub(keep_mydata);
y = data_meg.md(keep_megdata_sub);
c3_color = [0.8500 0.3250 0.0980];
plot(x, y, 'o', 'MarkerEdgeColor', c3_color, 'MarkerFaceColor', c3_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xylim_lo xylim_hi];
xax.TickValues = [xylim_lo (xylim_lo+xylim_hi)/2 xylim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [xylim_lo xylim_hi];
yax.TickValues = [xylim_lo (xylim_lo+xylim_hi)/2 xylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(xylim_lo, '%1.1f'), '', num2str(xylim_hi, '%1.1f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

legend([{'Equality'}, {'CA1'}, {'CA23dg'}, {'SUB'}], 'Location', 'northwest')
legend box off

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = {'Meg''s Data'; '(MD x 1e-3)'};
a.XLabel.String = {'Sophia''s Data'; '(MD x 1e-3)'};
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_equality_md_meg_bl_mrtrix3act'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_equality_md_meg_bl_mrtrix3act'), '-depsc')

hold off;



