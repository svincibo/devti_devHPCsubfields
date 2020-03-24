
clear all; close all; clc
format shortG

yc_color = [0 0.4470 0.7410];
oc_color = [0.4660 0.6740 0.1880];
a_color = [0.6350 0.0780 0.1840];

blprojectid = '5e5672430f7fa65e1d3c9621';

remove_outliers = 'no';
include = 'all'; % all, childrenonly
% outlier = [];

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/DevTI/devti_devHPCsubfields/';
% addpath(genpath(fullfile(rootDir, 'proj-5e5672430f7fa65e1d3c9621')));

% Read in behavioral data.
load(fullfile(rootDir, 'data.mat'))
beh_data_in_tbl = array2table(data, 'VariableNames', {'subID', 'cov_age', 'iq', 'gp_age', 'a', 'b', 'c', 'd', 'e'});

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, ['proj-' blprojectid]));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for t = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh_data_in_tbl.subID == str2num(grp_contents(t).name(6:7)))))
        
        % Display current sub ID.
        disp(grp_contents(t).name)
        
        % Update subject counter for when not all subjects are used/needed.
        sub_count = sub_count + 1;
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(t).folder, grp_contents(t).name, '/dt-raw.tag-snr*/*product.json'));
        
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        % Get SNR for this subject.
        data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
        
        for g = 1:size(data_snr_temp.SNRInB0_X_Y_Z, 1)
            
            get_temp(g) = str2num(data_snr_temp.SNRInB0_X_Y_Z{g});
            
        end
        
        % Get SNR.
        snr(sub_count) = min(get_temp);
        
        % Get subID.
        subID(sub_count) = str2num(grp_contents(t).name(5:7));
        
        % Get age group.
        group(sub_count) = beh_data_in_tbl.gp_age(find((beh_data_in_tbl.subID == str2num(grp_contents(t).name(6:7)))));
        
        %     % Get age in months.
        age(sub_count) = beh_data_in_tbl.cov_age(find((beh_data_in_tbl.subID == str2num(grp_contents(t).name(5:7)))));
        
        clear data_snr_temp get_temp
        
    end % end if exist
    
end % end t

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    subID = beh_data_in_tbl.subID;
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    snr = snr(~idx_outlier);
    group = group(~idx_outlier);
    age = age(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'SNR outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'SNR outlier retained.';
    
end

% Write out table for anova.
t_out = array2table(cat(2, subID', group', age', snr'), 'VariableNames', {'subID', 'group_age', 'cov_age', 'snr'});
writetable(t_out, fullfile(rootDir, 'supportFiles', ['devti_data_snr_' include '.csv']));

disp('Check for group differences in SNR.')
[~, tableout, ~] = anova1(snr', group', 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

disp('Check for correlation between age and SNR.')
[rho, p] = corr([snr' group']);
disp(['r = ' num2str(rho(1, 2)) ', p = ' num2str(p(1, 2)) '.'])

% Visualize: correlation.
figure(1)
plot(age(group==1)', snr(group==1)', 'LineStyle', 'none', 'MarkerEdgeColor', yc_color, 'MarkerFaceColor', yc_color, 'Marker', 'o', 'MarkerSize', 10) % use all subjects for this
hold on
% r_yc = plotcorr(age(group==1)', snr(group==1)', [], [], [], 'c');
plot(age(group==2)', snr(group==2)', 'LineStyle', 'none', 'MarkerEdgeColor', oc_color, 'MarkerFaceColor', oc_color, 'Marker', 'o', 'MarkerSize', 10)
% r_oc = plotcorr(age(group==2)', snr(group==2)', [], [], [], 'b');
[rc, pc] = plotcorr(age(group~=3)', snr(group~=3)', [], [], [], 'k');
if strcmp(include, 'all')
    plot(age(group==3)', snr(group==3)', 'r.', 'MarkerSize', 15)
    r_a = plotcorr(age(group==3)', snr(group==3)', [], [], [], 'r');
    r_all = plotcorr(age', snr', [], [], [], 'k');
    xlim_lo = min(age)-1;
    xlim_hi = max(age)+1;
    ylim_lo = min(snr)-1;
    ylim_hi = max(snr)+1;
else
    legend({'Children', 'Adolescents', ['Children, r = ' num2str(rc, '%.3f') 'p = ' num2str(pc, '%.3f')]}, 'Location', 'southwest')
    xlim_lo = min(age(group~=3))-1;
    xlim_hi = max(age(group~=3))+1;
    ylim_lo = min(snr(group~=3))-1;
    ylim_hi = max(snr(group~=3))+1;
end
legend('boxoff');
legend('off');

xlabel('Age (years)')
title('Correlations between Age and SNR')

capsize = 0;
marker = 'o';
markeredgecolor_h = [0 .73 .73]'; markeredgecolor_v = [.146 0 0]';
markerfacecolor_h = [0 .73 .73]; markerfacecolor_v = [.146 0 0];
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
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {num2str(xlim_lo, '%1.2f'), '', num2str(xlim_hi, '%1.2f')};
%     xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
% ylim_lo = 5; ylim_hi = 18;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_scatter_snr_' include]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_scatter_snr_' include]), '-depsc')

hold off;

% Visualize: group differences
figure(2)
hold on;
b1 = bar(1, nanmean(snr(group == 1)), 'FaceColor', yc_color, 'EdgeColor', yc_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(group == 1)) - nanstd(snr(group == 1)) nanmean(snr(group == 1)) + nanstd(snr(group == 1))], 'Color', yc_color)
b2 = bar(2, nanmean(snr(group == 2)), 'FaceColor', oc_color, 'EdgeColor', oc_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(group == 2)) - nanstd(snr(group == 2)) nanmean(snr(group == 2)) + nanstd(snr(group == 2))], 'Color', oc_color)
b3 = bar(3, nanmean(snr(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([3 3], [nanmean(snr(group == 3)) - nanstd(snr(group == 3)) nanmean(snr(group == 3)) + nanstd(snr(group == 3))], 'Color', a_color)

xlim_lo = min(age)-1;
xlim_hi = max(age)+1;
ylim_lo = min(snr)-1;
ylim_hi = max(snr)+1;
    
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

a.YLabel.String = 'SNR';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', ['plot_barplot_snr_' include]), '-dpng')
print(fullfile(rootDir, 'plots', 'eps', ['plot_barplot_snr_' include]), '-depsc')

hold off;
