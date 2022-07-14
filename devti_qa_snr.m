% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

remove_outliers = 'yes';
% Identify outliers to be removed.
outlier = [11 90];

include = 'all';

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Read in behavioral data.
load(fullfile(rootDir, 'supportFiles/data.mat'));
% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)

    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((data(:, 1)  == str2num(grp_contents(s).name(5:7)))))

        % Display current sub ID.
        disp(grp_contents(s).name)

        % Update subject counter for when not all subjects are used/needed.
        sub_count = sub_count + 1;

        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');

        % Get SNR for this subject.
        data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));

        % Get SNR in b0 images.
        b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});

        % Get mean SNR in X, Y, and Z directions.
        m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);

        % Get subID.
        subID(sub_count) = str2num(grp_contents(s).name(5:7));

        % Get training group.
        group(sub_count) = data(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))), 4);

        % Get age
        age(sub_count) = data(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))), 2);

        clear data_snr_temp get_temp

    end % end if exist

end % end t

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')

    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);

    % Remove outliers.
    subID = subID(~idx_outlier);
    b0 = b0(~idx_outlier);
    m = m(~idx_outlier);
    group = group(~idx_outlier);
    age = age(~idx_outlier);

    % Set figure note.
    ttlstr = 'SNR outlier removed.';

else

    % Set figure note.
    ttlstr = 'SNR outlier retained.';

end

% Write out table for anova.
t_out = array2table(cat(2, subID', group', m', b0'), 'VariableNames', {'subID', 'group', 'm', 'b0'});
writetable(t_out, fullfile(rootDir, 'supportFiles', 'devti_data_snr_indetail.csv'));

% Group differences test
disp('Are there b0 SNR differences among groups?')
[p, tableout, stats] = anova1(b0, group, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between children and adolescents?')
[h, p, ci stats] = ttest2(b0(group == 1), b0(group == 2));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between adolescents and adults?')
[h, p, ci stats] = ttest2(b0(group == 2), b0(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between children and adults?')
[h, p, ci stats] = ttest2(b0(group == 1), b0(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

disp('Are there weighted SNR differences among groups?')
[p, tableout, stats] = anova1(m, group, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between children and adolescents?')
[h, p, ci stats] = ttest2(m(group == 1), m(group == 2));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between adolescents and adults?')
[h, p, ci stats] = ttest2(m(group == 2), m(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between children and adults?')
[h, p, ci stats] = ttest2(m(group == 1), m(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Group differences plot: b0
snr = b0;
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
xtickvalues = [1 2 3];
alphablend = .8;
ylim_lo = 0;
ylim_hi = 40;

yc_color  = [50 180 100]/255;
oc_color = [50 100 180]/255;
a_color = [100 50 180]/255;

% Controls
b1 = bar(1, nanmean(snr(group == 1)), 'FaceColor', yc_color, 'EdgeColor', yc_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(group == 1)) - nanstd(snr(group == 1)) nanmean(snr(group == 1)) + nanstd(snr(group == 1))], 'Color', yc_color)
% Beginners
b2 = bar(2, nanmean(snr(group == 2)), 'FaceColor', oc_color, 'EdgeColor', oc_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(group == 2)) - nanstd(snr(group == 2)) nanmean(snr(group == 2)) + nanstd(snr(group == 2))], 'Color', oc_color)
% Experts
b3 = bar(3, nanmean(snr(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([3 3], [nanmean(snr(group == 3)) - nanstd(snr(group == 3)) nanmean(snr(group == 3)) + nanstd(snr(group == 3))], 'Color', a_color)

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

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR, b0 volumes';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_barplot_snr_b0_bygroup'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_barplot_snr_b0_bygroup'), '-depsc')

hold off;
clear snr

% Group differences plot: weighted
snr = m;
figure(2)
hold on;

% Controls
b1 = bar(1, nanmean(snr(group == 1)), 'FaceColor', yc_color, 'EdgeColor', yc_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(group == 1)) - nanstd(snr(group == 1)) nanmean(snr(group == 1)) + nanstd(snr(group == 1))], 'Color', yc_color)
% Beginners
b2 = bar(2, nanmean(snr(group == 2)), 'FaceColor', oc_color, 'EdgeColor', oc_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(group == 2)) - nanstd(snr(group == 2)) nanmean(snr(group == 2)) + nanstd(snr(group == 2))], 'Color', oc_color)
% Experts
b3 = bar(3, nanmean(snr(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([3 3], [nanmean(snr(group == 3)) - nanstd(snr(group == 3)) nanmean(snr(group == 3)) + nanstd(snr(group == 3))], 'Color', a_color)

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

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR, weighted volumes';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_barplot_snr_weighted_bygroup'), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', 'plot_barplot_snr_weighted_bygroup'), '-depsc')

hold off;

[r, p] = corrcoef(age, b0)
[r, p] = corrcoef(age, m)

