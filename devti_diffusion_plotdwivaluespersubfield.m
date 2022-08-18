% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

% Note bl-repository.
blprojectid = '5e5672430f7fa65e1d3c9621';

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Select measure.
wm = 'fa'; %'ad', 'rd', 'md'

% Load data.
load(fullfile(rootdir, 'supportFiles', ['devti_data_' wm '_20220705.mat']));

% Scale md values for analysis and visualization, if this is md.
if strcmp(wm, 'md')
    scalefactor = 1000;
elseif strcmp(wm, 'fa')
    scalefactor = 1;
end
m(:, 7:end) = array2table(table2array(m(:, 7:end)).*scalefactor);

% % Sort according to age. 
% age_idx = find(strcmp(m.Properties.VariableNames, 'age_cov') == 1);
% m = sortrows(m, age_idx);

% Select rois.
rois= m.Properties.VariableNames(7:end);

% Specify values for plotting figures.
capsize = 0;
marker = 'o';
linewidth = 4;
linestyle = 'none';
markersize = 10;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 24;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.02;

% Specify colors.
ca1 = [0 127 255]/255; % dark blue
ca23 = [240 128 128]/255; % light coral
dg = [147 112 219]/255; % medium purple
sub = [60 179 113]/255; %[64 224 208]/255; % turquoise

head = [204 204 0]/255; % light burnt yellow
body = [204 190 0]/255; % burnt yellow
tail =  [210 43 43]/255; % cadmium red

hip = [2 129 129]/255; % dark turquoise

yc_color  = [50 180 100]/255;
oc_color = [50 100 180]/255;
a_color = [100 50 180]/255;

for r = 1:length(rois)

    % Find the column that contains the data for the current roi and
    % extract the data.
    idx_roi = find(strcmp(m.Properties.VariableNames, rois{r}) == 1);
    roidata = table2array(m(:, idx_roi));

    figure(r); hold on;
    set(gcf,'Visible','on');
    set(0,'DefaultFigureVisible','on');

    % Plot data with age group color coded.
    s = gscatter(1:length(m.subID), roidata, m.age_group, [yc_color; oc_color; a_color]);
    s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
    s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
    s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;

    a = gca;
    if strcmp(wm, 'md')
        a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
    elseif strcmp(wm, 'fa')
        a.YLabel.String = {'Fractional Anisotropy (FA)'};
    end
    a.YLabel.FontSize = fontsize;
    a.YLabel.FontAngle = fontangle;

    a.XLabel.String = {'Subject'};
    a.XLabel.FontSize = fontsize;
    title(rois{r})

    % xaxis
    xlim_lo = 0; xlim_hi = 60;
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
    xax.TickDirection = 'out';
    xax.TickLength = [yticklength yticklength];
    xax.TickLabels = {num2str(xlim_lo, '%1.0f'), num2str((xlim_lo+xlim_hi)/2, '%1.0f'), num2str(xlim_hi, '%1.0f')};
    xax.FontName = fontname;
    xax.FontSize = fontsize;

    % yaxis
    if strcmp(wm, 'md')
        ylim_lo = 0.6; ylim_hi = 1.3;
    elseif strcmp(wm, 'fa')
        ylim_lo = 0; ylim_hi = 0.5;
    end
    yax = get(gca,'yaxis');
    yax.Limits = [ylim_lo ylim_hi];
    yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
    yax.TickDirection = 'out';
    yax.TickLength = [xticklength xticklength];
    yax.TickLabels = {num2str(ylim_lo, '%1.1f'), num2str((ylim_lo+ylim_hi)/2, '%1.1f'), num2str(ylim_hi, '%1.1f')};
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    yax.FontAngle = fontangle;

    % Add a gray band to indicate value range from prior studies.
    if strcmp(wm, 'md')
        %lowval = 0.9; highval = 1.5; % Langnes et al
        %lowval = 0.7; highval = 0.8; % Muller et al
        %lowval = 0.7; highval = 1.2; % Anblagan et al
        lowval = 0.7; highval = 1.2; %low/high for Muller and Anblagan
    elseif strcmp(wm, 'fa')
        %lowval = 0.2; highval = 0.3; % Muller et al
        %lowval = 0.08; highval = 0.15; % Anblagan et al
        lowval = 0.08; highval = 0.3; %low/high for Muller and Anblagan
    end
    patch([xlim_lo xlim_hi xlim_hi xlim_lo], [lowval lowval highval highval], [200 200 200]/255, 'FaceAlpha',0.2, 'EdgeColor','none')

    % legend off;
    legend({'Children', 'Adolescents', 'Adults'})
    box off;
    legend('box', 'off');
    legend('location', 'northeast');
    pbaspect([1 1 1])

    n = size(m, 1);
    print(fullfile(rootdir, 'plots', [wm '_' rois{r} '_n=' num2str(n)]), '-dpng')
    print(fullfile(rootdir, 'plots', 'eps', [wm '_' rois{r} '_n=' num2str(n)]), '-depsc')

end
