clear all; close all; clc
format shortG

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

blprojectid = '5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields/';
% addpath(genpath(fullfile(rootDir, 'proj-5e5672430f7fa65e1d3c9621')));

remove_outliers = 'yes';
include = 'all'; % all, childrenonly
outlier = [90];

%%%%%%%%%%%%%%% TESTING AREA %%%%%%%%%%%%%%%%

% % Use only the subjects that Meg used.
% data_meg = readtable(fullfile(rootDir, 'supportFiles/dti_subfields_long.csv'));
% sub_include = unique(data_meg.subnr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in behavioral data.
load(fullfile(rootDir, 'supportFiles/data.mat'))
data_beh = array2table(data, 'VariableNames', {'subID', 'cov_age', 'iq', 'gp_age', 'a', 'b', 'c', 'd', 'e'});

% Read in additional behavioral data, sex.
data_sex = readtable(fullfile(rootDir, 'supportFiles/devti_sex_all.csv'));

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, ['proj-' blprojectid]));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% One subject at a time.
sub_count = 0;
for t = 1:size(grp_contents, 1)
           
    % Only collect values for subjects that have both MRI and behavioral/demographic data.
    if ~isempty(find((data_beh.subID == str2num(grp_contents(t).name(6:7))))) 
    
    % Display current sub ID.
    disp(grp_contents(t).name)
    
    % Update subject counter for when not all subjects are used/needed.
    sub_count = sub_count + 1;
    
    % Get contents of the directory where the SNR values for this subject are stored.
    sub_contents_motion = dir(fullfile(grp_contents(t).folder, grp_contents(t).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
    % Remove the '.' and '..' files.
    sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
    
    % Get SNR for this subject.
    load([sub_contents_motion.folder filesep sub_contents_motion.name]);
    
    % Select only the translation/rotation parameters.
    mot = vertcat(xform(:).ecParams);
    mot = mot(:, 1:6); % xyz translation xyz rotation
    
    % Motion parameters represent the translation/rotation that must occur
    % to return the image at timepoint tn to the place that it was at timepoint
    % t0. Thus, to calculate instantaneouse parameters, we need a moving
    % difference.
    for m = 1:size(mot, 2)
        
        % Get moving difference for each parameter. Append row of zeros for t1, per convention (Power et al., 2014).
        % This step creates the fd timecourses for each motion parameter for each run in the motion struct and represents instantaneous movement.
        movingdifference(:, m) = [0 ; diff(mot(:, m), 1, 1)]';
        
    end
    
    % Get an overall fd for all 6 parameters for each run.
    % This step creates the fd summary statistic for all 6 parameters for each timepoint in a run for each run (e.g., scalar FD timecourse).
    motion(sub_count, :) = sum(abs(movingdifference), 2)';
    
    % Get subID.
    subID(sub_count) = str2num(grp_contents(t).name(5:7));
    
    % Get age group.
    group(sub_count) = data_beh.gp_age(find((data_beh.subID == str2num(grp_contents(t).name(6:7)))));
    
    %     % Get age in months.
    age(sub_count) = data_beh.cov_age(find((data_beh.subID == str2num(grp_contents(t).name(5:7)))));
    
    clear data_snr_temp get_temp
    
    end % end if exist
    
end % end t

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    motion = motion(~idx_outlier, :);
    group = group(~idx_outlier);
    age = age(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'Motion outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'Motion outlier retained.';
    
end

meanmotion = mean(motion, 2);

% Write out table for anova.
t_out = array2table(cat(2, subID', group', age', meanmotion), 'VariableNames', {'subID', 'group_age', 'cov_age', 'fd'});
writetable(t_out, fullfile(rootDir, 'supportFiles', ['devti_data_motion_' include '.csv']));

disp('Check for group differences in FD.')
[~, tableout, ~] = anova1(meanmotion, group', 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

disp('Check for correlation between age and FD.')
[rho, p] = corr([meanmotion group']);
disp(['r = ' num2str(rho(1, 2)) ', p = ' num2str(p(1, 2)) '.'])

%Visualize: correlation.
figure(1)
plot(age(group==1)', meanmotion(group==1)', 'LineStyle', 'none', 'MarkerEdgeColor', yc_color, 'MarkerFaceColor', yc_color, 'Marker', 'o', 'MarkerSize', 10) % use all subjects for this
hold on
plot(age(group==2)', meanmotion(group==2)', 'LineStyle', 'none', 'MarkerEdgeColor', oc_color, 'MarkerFaceColor', oc_color, 'Marker', 'o', 'MarkerSize', 10)
% [rc, pc] = plotcorr(age(group~=3)', meanmotion(group~=3), [], [], [], 'k:');
if strcmp(include, 'all')
    plot(age(group==3)', meanmotion(group==3), 'LineStyle', 'none', 'MarkerEdgeColor', a_color, 'MarkerFaceColor', a_color, 'Marker', 'o', 'MarkerSize', 10)
%     r_a = plotcorr(age(group==3)', meanmotion(group==3), [], [], [], 'k:');
    r_all = plotcorr(age', meanmotion, [], [], [], 'k');
    xlim_lo = min(age)-1;
    xlim_hi = max(age)+1;
    ylim_lo = min(meanmotion)-1;
    ylim_hi = max(meanmotion)+1;
else
    legend({'Children', 'Adolescents', ['Children, r = ' num2str(rc, '%.3f') 'p = ' num2str(pc, '%.3f')]}, 'Location', 'southwest')
    xlim_lo = min(age)-1;
    xlim_hi = max(age)+1;
    ylim_lo = min(meanmotion(group~=3))-1;
    ylim_hi = max(meanmotion(group~=3))+1;
end
legend('boxoff');
legend('off');

xlabel('Age (years)')
% title('Correlations between Age and FD')

capsize = 0;
marker = 'o';
% markeredgecolor_h = [0 .73 .73]'; markeredgecolor_v = [.146 0 0]';
% markerfacecolor_h = [0 .73 .73]; markerfacecolor_v = [.146 0 0];
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
xtickvalues = [1 2 3];
alphablend = .8;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 30];
xax.TickValues = [0 5 10 15 20 25 30];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'0', '5', '10', '15', '20', '25', '30'};
%     xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'Mean Framewise Displacement (FD)';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_scatter_motion_' include]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_scatter_motion_' include]), '-depsc')

hold off;

% Visualize: group differences
figure(2)
hold on;
b1 = bar(1, nanmean(meanmotion(group == 1)), 'FaceColor', yc_color, 'EdgeColor', yc_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(meanmotion(group == 1)) - nanstd(meanmotion(group == 1)) nanmean(meanmotion(group == 1)) + nanstd(meanmotion(group == 1))], 'Color', yc_color)
b2 = bar(2, nanmean(meanmotion(group == 2)), 'FaceColor', oc_color, 'EdgeColor', oc_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(meanmotion(group == 2)) - nanstd(meanmotion(group == 2)) nanmean(meanmotion(group == 2)) + nanstd(meanmotion(group == 2))], 'Color', oc_color)
b3 = bar(3, nanmean(meanmotion(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([3 3], [nanmean(meanmotion(group == 3)) - nanstd(meanmotion(group == 3)) nanmean(meanmotion(group == 3)) + nanstd(meanmotion(group == 3))], 'Color', a_color)

xlim_lo = min(age)-1;
xlim_hi = max(age)+1;
% ylim_lo = min(meanmotion)-1;
% ylim_hi = max(meanmotion)+1;
    
% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 3.5];
xax.TickValues = [1 2 3];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Children', 'Adolescents', 'Adults'};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'Mean Framewise Displacement (FD)';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_barplot_motion_all'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_barplot_motion_all'), '-depsc')

hold off;

% Manually record outliers. Include observations with unusually high FD and any above 2mm. 
% (0 indicates no outliers)
outliers.motion = 90;
% outliers.motion = subID(meanmotion>2);

save('devti_remove_motionoutliers.mat', 'outliers')
