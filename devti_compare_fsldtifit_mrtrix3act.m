% Age - Microstructural Relationships Among HPC Subfields

% For fits from FSLDTIFIT, use: devti_data_md.mat
% For fits from mrtrix3 act, use: devti_data_md_mrtrix3act.mat %NOTE: need to scale m by 1e3 on line 21 to match.

clear all; close all; clc
format long g

wm = {'fa', 'md'};

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields/';

% Bring in Meg's data.
data_meg = readtable(fullfile(rootDir, 'supportFiles/dti_subfields_long.csv'));
sub_include = unique(data_meg.subnr);

for w = 1:length(wm)
    
    if strcmp(wm{w}, 'fa')
        scale = 1;
    elseif strcmp(wm{w}, 'md')
        scale = 1000;
    end
    
    % Bring in my data: mrtrix3 act.
    load(fullfile(rootDir, ['supportFiles/devti_data_' wm{w} '_mrtrix3act.mat']))
    
    % Convert data to table for easier model specification.
    data = array2table(cat(2, transpose(sub), transpose(age), m), 'VariableNames', {'subID', 'age', roi{1, :}});
    
    % Find indices for subjects in my dataset that are also in Meg's dataset.
    keep_mydata = find(ismember(data.subID, data_meg.subnr));
    
    % Clearly set x and y values for ease.
    ca1_mrtrix3act = data.b_ca1(keep_mydata)*scale;
    ca23_mrtrix3act = data.b_ca23(keep_mydata)*scale;
    sub_mrtrix3act = data.b_sub(keep_mydata)*scale;
    
    clear data sub age sex
    
    % Bring in my data: FSLDTIFIT.
    load(fullfile(rootDir, ['supportFiles/devti_data_' wm{w} '.mat']))
    
    % Convert data to table for easier model specification.
    data = array2table(cat(2, transpose(sub), transpose(age), m), 'VariableNames', {'subID', 'age', roi{1, :}});
    
    % Find indices for subjects in my dataset that are also in Meg's dataset.
    keep_mydata = find(ismember(data.subID, data_meg.subnr));
    
    % Clearly set x and y values for ease.
    ca1_fsldtifit = data.b_ca1(keep_mydata);
    ca23_fsldtifit = data.b_ca23(keep_mydata);
    sub_fsldtifit = data.b_sub(keep_mydata);
    
    clear data sub age sex
    
    %% Visualize.
    
    figure(w)
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
    if strcmp(wm{w}, 'fa')
        xylim_lo = 0;
        xylim_hi = 0.5;
    elseif strcmp(wm{w}, 'md')
        xylim_lo = 0.5;
        xylim_hi = 1.4;
    end
    xtickvalues = [xylim_lo (xylim_lo+xylim_hi)/2 xylim_hi];
    
    plot(linspace(xylim_lo, xylim_hi, 2), linspace(xylim_lo, xylim_hi, 2), ':k')
    hold on;
    
    % Get md diff for ca1.
    x = ca1_mrtrix3act;
    y = ca1_fsldtifit;
    c1_color = [0 0.4470 0.7410];
    plot(x, y, 'o', 'MarkerEdgeColor', c1_color, 'MarkerFaceColor', c1_color)
    
    % Get md difffor ca23.
    x = ca23_mrtrix3act;
    y = ca23_fsldtifit;
    c2_color = [0.4940 0.1840 0.5560];
    plot(x, y, 'o', 'MarkerEdgeColor', c2_color, 'MarkerFaceColor', c2_color)
    
    % Get md difffor sub.
    x = sub_mrtrix3act;
    y = sub_fsldtifit;
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
    
    a.YLabel.String = {[ wm{w} ', FSLDTIFIT']};
    a.XLabel.String = {[wm{w} ', mrtrix3 act']};
    a.YLabel.FontSize = fontsize;
    pbaspect([1 1 1])
    
    print(fullfile(rootDir, 'plots', ['plot_equality_' wm{w} '_bl_fsldtifit_mrtrix3act']), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_equality_' wm{w} '_bl_fsldtifit_mrtrix3act']), '-depsc')
    
    hold off;
    
end



