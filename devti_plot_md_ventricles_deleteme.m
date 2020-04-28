% Age - Microstructural Relationships Among HPC Subfields

% This is a quick script written to plot the MD values in the ventricles of
% the DevTI_devHPCsubfields data.

clear all; close all; clc
format long g

wm = 'md';
proc = 'fsl'; %act, fsl

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields/';

% Get contents of the directory where the tract measures for this subject are stored.
sub_contents_rois = dir(fullfile(rootDir, blprojectid, ['sub-001/dt-raw.tag-tensor_metrics.' proc '*/' wm '*.nii.gz']));

% Remove the '.' and '..' files.
sub_contents_rois = sub_contents_rois(arrayfun(@(x) x.name(1), sub_contents_rois) ~= '.');

if strcmp(proc, 'act')
    scale = 1000;
elseif strcmp(proc, 'fsl')
    scale = 1;
end

for j = 1:size(sub_contents_rois)
    
    % Read in data for this subject and this tract.
    data_temp = niftiread(fullfile(sub_contents_rois(j).folder, sub_contents_rois(j).name));
    
    % Convert all zeros to NaN.
    data_temp(data_temp == 0) = NaN;
    
    % Get the indices of voxels that are a part of this ROI.
    idx = ~isnan(data_temp);
    
    % Select the MD values of voxels that are a part of this ROI.
    md{j} = data_temp(idx);
    
    % Assign group label that can be used to discriminate among MD values from different ROIs.
    roi(j) = j;
    
    clear data_temp
    
end % end j

%% Visualize.

figure
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

col = [0 0.4470 0.7410; 0.4940 0.1840 0.5560; 0.8500 0.3250 0.0980; 0.6 0.6 0; 0.4 0.6980 1; 1 0.6 0.2; 0.8 0.1 0.3];

hold on;

y_out = [];
for i = 1:size(md, 2)
    
    x = linspace(length(y_out)+1, length(y_out)+1+length(md{i}), length(md{i}));
    y = md{i}*scale;
    
    plot(x, y, 'o', 'MarkerEdgeColor', col(i, :), 'MarkerFaceColor', col(i, :))
    
    y_out = cat(1, y_out, y);
    
end

xlim_lo = 0; xlim_hi = length(y_out) + 1;
plot(linspace(xlim_lo, xlim_hi, 2), linspace(3, 3, 2), 'k')

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

ylim_lo = 0; ylim_hi = 4;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 3 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.1f'), '', '3.0', num2str(ylim_hi, '%1.1f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

legend([{'ROI1'}, {'ROI2'}, {'ROI3'}, {'ROI4'}], 'Location', 'northwest')
legend box off

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = {[wm ', ' proc]};
a.XLabel.String = {'voxel'};
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_scatter_' wm '_bl_ventricles_' proc]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_scatter_' wm '_bl_ventricles_' proc]), '-depsc')

hold off;




