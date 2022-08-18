% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

f = figure('visible', 'on');
set(gcf,'Visible','on');              
set(0,'DefaultFigureVisible','on');

accuracycutoff = 'no'; % selecting "no" going forward (7/8/2022) because includeBoolPerf and includeBookRT 
% exclude subjects who had poor performance due to not understanding the task, they had a complex way of 
% determining whether poor performance was a lack of understanding the task or just poor performance

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Load behavioral data.
load(fullfile(rootdir, 'supportFiles', 'devti_data_beh_20220705.mat'));
iq = readtable(fullfile(rootdir, 'supportFiles', 'devti_iq_all.xlsx')); 
iq.Properties.VariableNames(1) = {'subID'}; iq.Properties.VariableNames(2) = {'iq'};

% Bring in sex data.
sex = struct2table(load(fullfile(rootdir, 'supportFiles', 'devti_sex_all.mat')));

% Rectify subIDs.
idx = find(ismember(iq.subID, beh.subID));
iq = iq(idx, :); clear idx;
idx = find(ismember(sex.subID, beh.subID));
sex = sex(idx, :); clear idx;

% COMBINE BEHAVIORAL AND DIFFUSION DATA
% % Load diffusion data.
% wm = 'fa';
% load(fullfile(rootdir, 'supportFiles', ['devti_data_' wm '_20220705.mat']));
% mri = m(:, [1 4:end]); clear m;
% % Scale md values for analysis and visualization.
% if strcmp(wm, 'md')
%     d(:, 5:end) = array2table(table2array(d(:, 5:end)).*1000);
% end
% 
% % Rectify subIDs.
% idx = find(ismember(mri.subID, beh.subID));
% mri = mri(idx, :); clear idx;
% idx = find(ismember(beh.subID, mri.subID));
% beh = beh(idx, :); clear idx;

% 
% % Concatenate beh and mri data into one table.
% d = [beh(:, 1:3) mri(:, 2:4) beh(:, 4:end) iq(:, 2) mri(:, 5:end)];

% USE ONLY BEHAVIORAL DATA
% Concatenate beh and mri data into one table.
d = [beh(:, 1:3) sex(:, 1) beh(:, 4:end) iq(:, 2)];

% Remove the random spaces and other issues in the column names.
d.Properties.VariableNames = strrep(d.Properties.VariableNames, ' ', '');
d.Properties.VariableNames = strrep(d.Properties.VariableNames, '|', '');
d.Properties.VariableNames = strrep(d.Properties.VariableNames, '&', '');

% Make things easier.
d.assoc = d.task1_acc_dirperf4;
%d.assoc = mean([d.task1_acc_dirperf1 d.task1_acc_dirperf2 d.task1_acc_dirperf3 d.task1_acc_dirperf4], 2);
d.infer = d.task1_acc_ACperf;

% Tell Matlab that sex and age group are categorical variables.
d.sex = categorical(d.sex);
% d.group = categorical(d.group); % not defined as categorical because used
% for indexing and not part of an anova performed within this code
clear beh; % mri;

% % Remove subjects that were below chance on any relational integration task,
% % i.e., rel-0, rel-1, or rel-2.
% idx1 = d.subID(find(d.assoc <= .50 | d.infer <= .50)); 
% idx4 = d.subID(find(d.task4_rel0acc <= .33 | d.task4_rel1acc <= .33 | d.task4_rel2acc <= .33)); %there were none
% idx5 = d.subID(find(d.task5_GenDIRacc <= .5 | d.task5_GenINDacc <= .50));

capsize = 0;
marker = 'o';
linewidth = 4;
linestyle = 'none';
markersize = 14;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 24;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.02;

%% Task 1 and Task 3: Final Direct and Inference, Schlichting et al., 2017
% Remove data identified as outliers by Schlichting et al., 2017. 
temp = load(fullfile(rootdir, 'taskperf', 'task1_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp removesubid;

% % % Remove data from subjects with accuracy < .05.
if strcmp(accuracycutoff, 'yes')
    d2.assoc(d2.assoc < .33) = NaN;
    d2.infer(d2.infer < .33) = NaN;
end

m1 = nanmean(d2.assoc(find(d2.group == 1))); n1 = length(d2.assoc(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.assoc(find(d2.group == 1)))/sqrt(n1));
m2 = nanmean(d2.infer(find(d2.group == 1))); n2 = length(d2.infer(find(d2.group == 1))); sem2 = 1.96*(nanstd(d2.infer(find(d2.group == 1)))/sqrt(n1));

m3 = nanmean(d2.assoc(find(d2.group == 2))); n3 = length(d2.assoc(find(d2.group == 2))); sem3 = 1.96*(nanstd(d2.assoc(find(d2.group == 2)))/sqrt(n2));
m4 = nanmean(d2.infer(find(d2.group == 2))); n4 = length(d2.infer(find(d2.group == 2))); sem4 = 1.96*(nanstd(d2.infer(find(d2.group == 2)))/sqrt(n2));

m5 = nanmean(d2.assoc(find(d2.group == 3))); n5 = length(d2.assoc(find(d2.group == 3))); sem5 = 1.96*(nanstd(d2.assoc(find(d2.group == 3)))/sqrt(n3));
m6 = nanmean(d2.infer(find(d2.group == 3))); n6 = length(d2.infer(find(d2.group == 3))); sem6 = 1.96*(nanstd(d2.infer(find(d2.group == 3)))/sqrt(n3));

figure(1); hold on;
b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([2], [m2]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([4], [m3]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([5], [m4]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([7], [m5]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([8], [m6]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;

p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
p = plot([4 4], [m3+sem3 m3-sem3]); p(1).Color = 'k';
p = plot([5 5], [m4+sem4 m4-sem4]); p(1).Color = 'k';
p = plot([7 7], [m5+sem5 m5-sem5]); p(1).Color = 'k';
p = plot([8 8], [m6+sem6 m6-sem6]); p(1).Color = 'k';

plot([0 9], [0.33 0.33], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 9];
xax.TickValues = [1.5 4.5 7.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n2)], ['adult, N=' num2str(n3)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 0.2 0.4 0.6 0.8 1];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str(0.2, '%1.2f'), num2str(0.4, '%1.2f'), num2str(0.6, '%1.2f'), num2str(0.8, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Accuracy (proportion correct)', 'FontAngle', 'italic')
   
% legend off;
legend({'Direct', 'Inference'})
box off; 
legend('box', 'off');
legend('location', 'eastoutside');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', 'devti_directinference_repfig4'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_directinference_repfig4'), '-depsc')

hold off;

[H, P, CI, stats] = ttest(d2.assoc(find(d2.group == 1)), d2.infer(find(d2.group == 1)));
[H, P, CI, stats] = ttest(d2.assoc(find(d2.group == 2)), d2.infer(find(d2.group == 2)));
[H, P, CI, stats] = ttest(d2.assoc(find(d2.group == 3)), d2.infer(find(d2.group == 3)));

[H, P, CI, stats] = ttest2(d2.assoc(find(d2.group == 1)), d2.assoc(find(d2.group == 2)));
[H, P, CI, stats] = ttest2(d2.infer(find(d2.group == 1)), d2.infer(find(d2.group == 2)));

% %% Task 2: Iowa Gambling Task
% % removesub = idx1;
% d2 = d;%(~ismember(d.subID, removesub), :);
% 
% m1 = nanmean(d2.task3_meanacc(find(d2.group == 1))); n1 = length(d2.task3_meanacc(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.task3_meanacc(find(d2.group == 1)))/sqrt(n1));
% m2 = nanmean(d2.task3_meanacc(find(d2.group == 2))); n2 = length(d2.task3_meanacc(find(d2.group == 2))); sem2 = 1.96*(nanstd(d2.task3_meanacc(find(d2.group == 2)))/sqrt(n2));
% m3 = nanmean(d2.task3_meanacc(find(d2.group == 3))); n3 = length(d2.task3_meanacc(find(d2.group == 3))); sem3 = 1.96*(nanstd(d2.task3_meanacc(find(d2.group == 3)))/sqrt(n3));
% 
% figure(2); hold on;
% b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([2], [m2]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([3], [m3]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% 
% p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
% p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
% p = plot([3 3], [m3+sem3 m3-sem3]); p(1).Color = 'k';
% 
% % plot([0 12], [0.50 0.50], 'k:');
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [0.5 3.5];
% xax.TickValues = [1 2 3];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n2)], ['adult, N=' num2str(n3)]};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 1;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [0 0.5 1];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
% 
% ylabel('Accuracy (proportion correct)', 'FontAngle', 'italic')
%    
% legend off;
% % legend({'DIR', 'IND'})
% box off; 
% % legend('box', 'off');
% % legend('location', 'eastoutside');
% pbaspect([1 1 1])
% 
% print(fullfile(rootdir, 'plots', 'devti_iowagamblingtask'), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', 'devti_iowagamblingtask'), '-depsc')
% 
% hold off;
% 
% %% Task 4: Relational Reasoning
% % Remove data identified as outliers by Schlichting et al., 2017. 
% temp = load(fullfile(rootdir, 'taskperf', 'task4_behavior_n=78.mat'));
% removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
% removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
% removesubid = unique([removetemp1 removetemp2]);
% clear temp removetemp1 removetemp2;
% d2 = d(~ismember(d.subID, removesubid), :); clear temp;
% 
% if strcmp(accuracycutoff, 'yes')
%     d2.task4_rel0acc(d2.task4_rel0acc < .33) = NaN;
%     d2.task4_rel1acc(d2.task4_rel1acc < .33) = NaN;
%     d2.task4_rel2acc(d2.task4_rel2acc < .33) = NaN;
% end
% 
% m1 = nanmean(d2.task4_rel0acc(find(d2.group == 1))); n1 = length(d2.task4_rel0acc(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.task4_rel0acc(find(d2.group == 1)))/sqrt(n1));
% m2 = nanmean(d2.task4_rel1acc(find(d2.group == 1))); n2 = length(d2.task4_rel1acc(find(d2.group == 1))); sem2 = 1.96*(nanstd(d2.task4_rel1acc(find(d2.group == 1)))/sqrt(n2));
% m3 = nanmean(d2.task4_rel2acc(find(d2.group == 1))); n3 = length(d2.task4_rel2acc(find(d2.group == 1))); sem3 = 1.96*(nanstd(d2.task4_rel2acc(find(d2.group == 1)))/sqrt(n3));
% 
% m4 = nanmean(d2.task4_rel0acc(find(d2.group == 2))); n4 = length(d2.task4_rel0acc(find(d2.group == 2))); sem4 = 1.96*(nanstd(d2.task4_rel0acc(find(d2.group == 2)))/sqrt(n4));
% m5 = nanmean(d2.task4_rel1acc(find(d2.group == 2))); n5 = length(d2.task4_rel1acc(find(d2.group == 2))); sem5 = 1.96*(nanstd(d2.task4_rel1acc(find(d2.group == 2)))/sqrt(n5));
% m6 = nanmean(d2.task4_rel2acc(find(d2.group == 2))); n6 = length(d2.task4_rel2acc(find(d2.group == 2))); sem6 = 1.96*(nanstd(d2.task4_rel2acc(find(d2.group == 2)))/sqrt(n6));
% 
% m7 = nanmean(d2.task4_rel0acc(find(d2.group == 3))); n7 = length(d2.task4_rel0acc(find(d2.group == 3))); sem7 = 1.96*(nanstd(d2.task4_rel0acc(find(d2.group == 3)))/sqrt(n7));
% m8 = nanmean(d2.task4_rel1acc(find(d2.group == 3))); n8 = length(d2.task4_rel1acc(find(d2.group == 3))); sem8 = 1.96*(nanstd(d2.task4_rel1acc(find(d2.group == 3)))/sqrt(n8));
% m9 = nanmean(d2.task4_rel2acc(find(d2.group == 3))); n9 = length(d2.task4_rel2acc(find(d2.group == 3))); sem9 = 1.96*(nanstd(d2.task4_rel2acc(find(d2.group == 3)))/sqrt(n9));
% 
% figure(3); hold on;
% b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([2], [m2]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.5;
% b = bar([3], [m3]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% b = bar([5], [m4]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([6], [m5]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.3;
% b = bar([7], [m6]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% b = bar([9], [m7]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([10], [m8]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.5;
% b = bar([11], [m9]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
% p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
% p = plot([3 3], [m3+sem3 m3-sem3]); p(1).Color = 'k';
% 
% p = plot([5 5], [m4+sem4 m4-sem4]); p(1).Color = 'k';
% p = plot([6 6], [m5+sem5 m5-sem5]); p(1).Color = 'k';
% p = plot([7 7], [m6+sem6 m6-sem6]); p(1).Color = 'k';
% 
% p = plot([9 9], [m7+sem7 m7-sem7]); p(1).Color = 'k';
% p = plot([10 10], [m8+sem8 m8-sem8]); p(1).Color = 'k';
% p = plot([11 11], [m9+sem9 m9-sem9]); p(1).Color = 'k';
% 
% plot([0 12], [0.33 0.33], 'k:');
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [0 12];
% xax.TickValues = [2 6 10];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n4)], ['adult, N=' num2str(n7)]};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0.7; ylim_hi = 1;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [0.7 0.845 1];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
% 
% ylabel('Accuracy (proportion correct)', 'FontAngle', 'italic')
%    
% % legend off;
% legend({'Rel 0', 'Rel 1', 'Rel 2'})
% box off; 
% legend('box', 'off');
% legend('location', 'eastoutside');
% pbaspect([1 1 1])
% 
% print(fullfile(rootdir, 'plots', 'devti_relationalreasoning_acc'), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalreasoning_acc'), '-depsc')
% 
% [H, P, CI, stats] = ttest2(d2.task4_rel0acc(find(d2.group == 1)), d2.task4_rel0acc(find(d2.group == 2)))
% [H, P, CI, stats] = ttest2(d2.task4_rel1acc(find(d2.group == 1)), d2.task4_rel1acc(find(d2.group == 2)))
% [H, P, CI, stats] = ttest2(d2.task4_rel2acc(find(d2.group == 1)), d2.task4_rel2acc(find(d2.group == 2)))
% [H, P, CI, stats] = ttest2(d2.task4_rel2acc(find(d2.group == 2)), d2.task4_rel2acc(find(d2.group == 3)))
% 
% hold off;

% %% Task 5: Relational Inference
% Remove data identified as outliers by Schlichting et al., 2017. 
temp = load(fullfile(rootdir, 'taskperf', 'task5_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp;

% % Remove data from subjects with accuracy < .5.
if strcmp(accuracycutoff, 'yes')
    d2.task5_GenDIRacc(d2.task5_GenDIRacc < .5) = NaN;
    d2.task5_GenINDacc(d2.task5_GenINDacc < .5) = NaN;
end

m1 = nanmean(d2.task5_GenDIRacc(find(d2.group == 1))); n1 = length(d2.task5_GenDIRacc(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.task5_GenDIRacc(find(d2.group == 1)))/sqrt(n1));
m2 = nanmean(d2.task5_GenINDacc(find(d2.group == 1))); n2 = length(d2.task5_GenINDacc(find(d2.group == 1))); sem2 = 1.96*(nanstd(d2.task5_GenINDacc(find(d2.group == 1)))/sqrt(n1));

m3 = nanmean(d2.task5_GenDIRacc(find(d2.group == 2))); n3 = length(d2.task5_GenDIRacc(find(d2.group == 2))); sem3 = 1.96*(nanstd(d2.task5_GenDIRacc(find(d2.group == 2)))/sqrt(n2));
m4 = nanmean(d2.task5_GenINDacc(find(d2.group == 2))); n4 = length(d2.task5_GenINDacc(find(d2.group == 2))); sem4 = 1.96*(nanstd(d2.task5_GenINDacc(find(d2.group == 2)))/sqrt(n2));

m5 = nanmean(d2.task5_GenDIRacc(find(d2.group == 3))); n5 = length(d2.task5_GenDIRacc(find(d2.group == 3))); sem5 = 1.96*(nanstd(d2.task5_GenDIRacc(find(d2.group == 3)))/sqrt(n3));
m6 = nanmean(d2.task5_GenINDacc(find(d2.group == 3))); n6 = length(d2.task5_GenINDacc(find(d2.group == 3))); sem6 = 1.96*(nanstd(d2.task5_GenINDacc(find(d2.group == 3)))/sqrt(n3));

figure(4); hold on;
b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([2], [m2]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([4], [m3]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([5], [m4]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([7], [m5]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([8], [m6]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;

p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
p = plot([4 4], [m3+sem3 m3-sem3]); p(1).Color = 'k';
p = plot([5 5], [m4+sem4 m4-sem4]); p(1).Color = 'k';
p = plot([7 7], [m5+sem5 m5-sem5]); p(1).Color = 'k';
p = plot([8 8], [m6+sem6 m6-sem6]); p(1).Color = 'k';

plot([0 9], [0.5 0.5], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 9];
xax.TickValues = [1.5 4.5 7.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n2)], ['adult, N=' num2str(n3)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 0.2 0.4 0.6 0.8 1];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str(0.2, '%1.2f'), num2str(0.4, '%1.2f'), num2str(0.6, '%1.2f'), num2str(0.8, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Accuracy (proportion correct)', 'FontAngle', 'italic')
   
% legend off;
legend({'Direct', 'Inference'})
box off; 
legend('box', 'off');
legend('location', 'eastoutside');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', 'devti_relationalinference_general_acc'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalinference_general_acc'), '-depsc')

hold off;

figure(5); hold on;
m1 = nanmean(d2.task5_GenDIRacc); n1 = length(d2.task5_GenDIRacc); sem1 = 1.96*(nanstd(d2.task5_GenDIRacc)/sqrt(n1));

m3 = nanmean(d2.task5_GenINDacc); n3 = length(d2.task5_GenINDacc); sem3 = 1.96*(nanstd(d2.task5_GenINDacc)/sqrt(n3));

b = bar([1.5], [m1]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([4.5], [m3]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;

p = plot([1.5 1.5], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([4.5 4.5], [m3+sem3 m3-sem3]); p(1).Color = 'k';

plot([0 9], [0.5 0.5], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 6];
xax.TickValues = [1.5 4.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['Gen Direct, N=' num2str(n1)], ['Gen Inference, N=' num2str(n3)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 0.2 0.4 0.6 0.8 1];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str(0.2, '%1.2f'), num2str(0.4, '%1.2f'), num2str(0.6, '%1.2f'), num2str(0.8, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Accuracy (proportion correct)', 'FontAngle', 'italic')
   
% legend off;
legend({'Direct', 'Inference'})
box off; 
legend('box', 'off');
legend('location', 'eastoutside');
pbaspect([2 1 2])

print(fullfile(rootdir, 'plots', 'devti_relationalinference_general_acc_maineffectoftask'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalinference_general_acc_maineffectoftask'), '-depsc')

hold off;

[H, P, CI, stats] = ttest(d2.task5_GenDIRacc(find(d2.group == 1)), d2.task5_GenINDacc(find(d2.group == 1)));
[H, P, CI, stats] = ttest(d2.task5_GenDIRacc(find(d2.group == 2)), d2.task5_GenINDacc(find(d2.group == 2)));
[H, P, CI, stats] = ttest(d2.task5_GenDIRacc(find(d2.group == 3)), d2.task5_GenINDacc(find(d2.group == 3)));

[H, P, CI, stats] = ttest2(d2.task5_GenINDacc(find(d2.group == 1)), d2.task5_GenINDacc(find(d2.group == 2)))


%% Task 1 and Task 3: Final Direct and Inference, Schlichting et al., 2017
d.assocrt = nanmean([d.task1_rt_dirrt1 d.task1_rt_dirrt2 d.task1_rt_dirrt3 d.task1_rt_dirrt4], 2);
d.inferrt = d.task1_rt_ACrt;

% % Remove data from subjects with accuracy < .05.
if strcmp(accuracycutoff, 'yes')
    d2.assocrt(d2.assoc < .5) = NaN;
    d2.inferrt(d2.infer < .5) = NaN;
end

% Remove data identified as outliers by Schlichting et al., 2017. 
temp = load(fullfile(rootdir, 'taskperf', 'task1_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp;

m1 = nanmean(d2.assocrt(find(d2.group == 1))); n1 = length(d2.assocrt(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.assocrt(find(d2.group == 1)))/sqrt(n1));
m2 = nanmean(d2.inferrt(find(d2.group == 1))); n2 = length(d2.inferrt(find(d2.group == 1))); sem2 = 1.96*(nanstd(d2.inferrt(find(d2.group == 1)))/sqrt(n1));

m3 = nanmean(d2.assocrt(find(d2.group == 2))); n3 = length(d2.assocrt(find(d2.group == 2))); sem3 = 1.96*(nanstd(d2.assocrt(find(d2.group == 2)))/sqrt(n2));
m4 = nanmean(d2.inferrt(find(d2.group == 2))); n4 = length(d2.inferrt(find(d2.group == 2))); sem4 = 1.96*(nanstd(d2.inferrt(find(d2.group == 2)))/sqrt(n2));

m5 = nanmean(d2.assocrt(find(d2.group == 3))); n5 = length(d2.assocrt(find(d2.group == 3))); sem5 = 1.96*(nanstd(d2.assocrt(find(d2.group == 3)))/sqrt(n3));
m6 = nanmean(d2.inferrt(find(d2.group == 3))); n6 = length(d2.inferrt(find(d2.group == 3))); sem6 = 1.96*(nanstd(d2.inferrt(find(d2.group == 3)))/sqrt(n3));

figure(6); hold on;
b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([2], [m2]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([4], [m3]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([5], [m4]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([7], [m5]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([8], [m6]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;

p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
p = plot([4 4], [m3+sem3 m3-sem3]); p(1).Color = 'k';
p = plot([5 5], [m4+sem4 m4-sem4]); p(1).Color = 'k';
p = plot([7 7], [m5+sem5 m5-sem5]); p(1).Color = 'k';
p = plot([8 8], [m6+sem6 m6-sem6]); p(1).Color = 'k';

% plot([0 9], [0.5 0.5], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 9];
xax.TickValues = [1.5 4.5 7.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n2)], ['adult, N=' num2str(n3)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 6000;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 2000 4000 6000];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str(2000, '%1.0f'), num2str(4000, '%1.0f'), num2str(6000, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
   
% legend off;
legend({'Direct', 'Inference'})
box off; 
legend('box', 'off');
legend('location', 'eastoutside');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', 'devti_directinference_repfig4_rt'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_directinference_repfig4_rt'), '-depsc')

hold off;

figure(7); hold on;
vec = [d2.assocrt(find(d2.group == 1)); d2.inferrt(find(d2.group == 1))];
m1 = nanmean(vec); n1 = length(vec)/2; sem1 = 1.96*(nanstd(vec)/sqrt(n1));

vec = [d2.assocrt(find(d2.group == 2)); d2.inferrt(find(d2.group == 2))];
m3 = nanmean(vec); n3 = length(vec)/2; sem3 = 1.96*(nanstd(vec)/sqrt(n3));

vec = [d2.assocrt(find(d2.group == 3)); d2.inferrt(find(d2.group == 3))];
m5 = nanmean(vec); n5 = length(vec)/2; sem5 = 1.96*(nanstd(vec)/sqrt(n5));

b = bar([1.5], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([4.5], [m3]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([7.5], [m5]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;

p = plot([1.5 1.5], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([4.5 4.5], [m3+sem3 m3-sem3]); p(1).Color = 'k';
p = plot([7.5 7.5], [m5+sem5 m5-sem5]); p(1).Color = 'k';

% plot([0 9], [0.5 0.5], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 9];
xax.TickValues = [1.5 4.5 7.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n3)], ['adult, N=' num2str(n5)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 6000;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 2000 4000 6000];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str(2000, '%1.0f'), num2str(4000, '%1.0f'), num2str(6000, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
   
legend off;
% legend({'Direct', 'Inference'})
box off; 
% legend('box', 'off');
% legend('location', 'eastoutside');
pbaspect([2 1 2])

print(fullfile(rootdir, 'plots', 'devti_directinference_repfig4_rt_maineffectgroup'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_directinference_repfig4_rt_t_maineffectgroup'), '-depsc')

hold off;

figure(8); hold on;
m1 = nanmean(d2.assocrt); n1 = length(d2.assocrt); sem1 = 1.96*(nanstd(d2.assocrt)/sqrt(n1));

m3 = nanmean(d2.inferrt); n3 = length(d2.inferrt); sem3 = 1.96*(nanstd(d2.inferrt)/sqrt(n3));

b = bar([1.5], [m1]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([4.5], [m3]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;

p = plot([1.5 1.5], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([4.5 4.5], [m3+sem3 m3-sem3]); p(1).Color = 'k';

% plot([0 9], [0.5 0.5], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 6];
xax.TickValues = [1.5 4.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['AB Direct, N=' num2str(n1)], ['AC Inference, N=' num2str(n3)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 6000;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 2000 4000 6000];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str(2000, '%1.0f'), num2str(4000, '%1.0f'), num2str(6000, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
   
% legend off;
legend({'Direct', 'Inference'})
box off; 
legend('box', 'off');
legend('location', 'eastoutside');
pbaspect([2 1 2])

print(fullfile(rootdir, 'plots', 'devti_directinference_repfig4_rt_maineffecttask'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_directinference_repfig4_rt_t_maineffecttask'), '-depsc')

hold off;

[H, P, CI, stats] = ttest(d2.assocrt(find(d2.group == 1)), d2.inferrt(find(d2.group == 1)))
[H, P, CI, stats] = ttest(d2.assocrt(find(d2.group == 2)), d2.inferrt(find(d2.group == 2)))
[H, P, CI, stats] = ttest(d2.assocrt(find(d2.group == 3)), d2.inferrt(find(d2.group == 3)))

% %% Task 4: Relational Reasoning
% % Remove data identified as outliers by Schlichting et al., 2017. 
% temp = load(fullfile(rootdir, 'taskperf', 'task4_behavior_n=78.mat'));
% removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
% removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
% removesubid = unique([removetemp1 removetemp2]);
% clear temp removetemp1 removetemp2;
% d2 = d(~ismember(d.subID, removesubid), :); clear temp;
% 
% if strcmp(accuracycutoff, 'yes')
%     d2.task4_rel0rt(d2.task4_rel0acc < .33) = NaN;
%     d2.task4_rel1rt(d2.task4_rel1acc < .33) = NaN;
%     d2.task4_rel2rt(d2.task4_rel2acc < .33) = NaN;
% end
% 
% m1 = nanmean(d2.task4_rel0rt(find(d2.group == 1))); n1 = length(d2.task4_rel0rt(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.task4_rel0rt(find(d2.group == 1)))/sqrt(n1));
% m2 = nanmean(d2.task4_rel1rt(find(d2.group == 1))); n2 = length(d2.task4_rel1rt(find(d2.group == 1))); sem2 = 1.96*(nanstd(d2.task4_rel1rt(find(d2.group == 1)))/sqrt(n2));
% m3 = nanmean(d2.task4_rel2rt(find(d2.group == 1))); n3 = length(d2.task4_rel2rt(find(d2.group == 1))); sem3 = 1.96*(nanstd(d2.task4_rel2rt(find(d2.group == 1)))/sqrt(n3));
% 
% m4 = nanmean(d2.task4_rel0rt(find(d2.group == 2))); n4 = length(d2.task4_rel0rt(find(d2.group == 2))); sem4 = 1.96*(nanstd(d2.task4_rel0rt(find(d2.group == 2)))/sqrt(n4));
% m5 = nanmean(d2.task4_rel1rt(find(d2.group == 2))); n5 = length(d2.task4_rel1rt(find(d2.group == 2))); sem5 = 1.96*(nanstd(d2.task4_rel1rt(find(d2.group == 2)))/sqrt(n5));
% m6 = nanmean(d2.task4_rel2rt(find(d2.group == 2))); n6 = length(d2.task4_rel2rt(find(d2.group == 2))); sem6 = 1.96*(nanstd(d2.task4_rel2rt(find(d2.group == 2)))/sqrt(n6));
% 
% m7 = nanmean(d2.task4_rel0rt(find(d2.group == 3))); n7 = length(d2.task4_rel0rt(find(d2.group == 3))); sem7 = 1.96*(nanstd(d2.task4_rel0rt(find(d2.group == 3)))/sqrt(n7));
% m8 = nanmean(d2.task4_rel1rt(find(d2.group == 3))); n8 = length(d2.task4_rel1rt(find(d2.group == 3))); sem8 = 1.96*(nanstd(d2.task4_rel1rt(find(d2.group == 3)))/sqrt(n8));
% m9 = nanmean(d2.task4_rel2rt(find(d2.group == 3))); n9 = length(d2.task4_rel2rt(find(d2.group == 3))); sem9 = 1.96*(nanstd(d2.task4_rel2rt(find(d2.group == 3)))/sqrt(n9));
% 
% figure(9); hold on;
% b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([2], [m2]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.5;
% b = bar([3], [m3]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% b = bar([5], [m4]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([6], [m5]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.3;
% b = bar([7], [m6]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% b = bar([9], [m7]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([10], [m8]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.5;
% b = bar([11], [m9]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
% p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
% p = plot([3 3], [m3+sem3 m3-sem3]); p(1).Color = 'k';
% 
% p = plot([5 5], [m4+sem4 m4-sem4]); p(1).Color = 'k';
% p = plot([6 6], [m5+sem5 m5-sem5]); p(1).Color = 'k';
% p = plot([7 7], [m6+sem6 m6-sem6]); p(1).Color = 'k';
% 
% p = plot([9 9], [m7+sem7 m7-sem7]); p(1).Color = 'k';
% p = plot([10 10], [m8+sem8 m8-sem8]); p(1).Color = 'k';
% p = plot([11 11], [m9+sem9 m9-sem9]); p(1).Color = 'k';
% 
% % plot([0 12], [0.33 0.33], 'k:');
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [0 12];
% xax.TickValues = [2 6 10];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n4)], ['adult, N=' num2str(n7)]};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 6000;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [0 3000 6000];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str((ylim_lo+ylim_hi)/2, '%1.0f'), num2str(ylim_hi, '%1.0f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
% 
% ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
%    
% % legend off;
% legend({'Rel 0', 'Rel 1', 'Rel 2'})
% box off; 
% legend('box', 'off');
% legend('location', 'eastoutside');
% pbaspect([1 1 1])
% 
% print(fullfile(rootdir, 'plots', 'devti_relationalreasoning_rt'), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalreasoning_rt'), '-depsc')
% 
% hold off;
% 
% figure(10); hold on;
% m1 = nanmean(d2.task4_rel0rt); n1 = length(d2.task4_rel0rt); sem1 = 1.96*(nanstd(d2.task4_rel0rt)/sqrt(n1));
% m2 = nanmean(d2.task4_rel1rt); n2 = length(d2.task4_rel1rt); sem2 = 1.96*(nanstd(d2.task4_rel1rt)/sqrt(n2));
% m3 = nanmean(d2.task4_rel2rt); n3 = length(d2.task4_rel2rt); sem3 = 1.96*(nanstd(d2.task4_rel2rt)/sqrt(n3));
% 
% figure(10); hold on;
% b = bar([1], [m1]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([2], [m2]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.5;
% b = bar([3], [m3]); color = [128 128 128]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
% 
% p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
% p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
% p = plot([3 3], [m3+sem3 m3-sem3]); p(1).Color = 'k';
% 
% % plot([0 12], [0.33 0.33], 'k:');
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [0 4];
% xax.TickValues = [1 2 3];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {['Rel-0, N=' num2str(n1)], ['Rel-1, N=' num2str(n2)], ['Rel-2, N=' num2str(n3)]};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 6000;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [0 3000 6000];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str((ylim_lo+ylim_hi)/2, '%1.0f'), num2str(ylim_hi, '%1.0f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
% 
% ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
%    
% % legend off;
% legend({'Rel-0', 'Rel-1', 'Rel-2'})
% box off; 
% legend('box', 'off');
% legend('location', 'eastoutside');
% pbaspect([2 1 2])
% 
% print(fullfile(rootdir, 'plots', 'devti_relationalreasoning_rt_maineffectoftask'), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalreasoning_rt_maineffectoftask'), '-depsc')
% 
% hold off;
% 
% [H, P, CI, stats] = ttest2(d2.task4_rel2rt(find(d2.group == 1)), d2.task4_rel2rt(find(d2.group == 2)))
% 
% figure(10); hold on;
% vec = [d2.task4_rel0rt(find(d2.group == 1)); d2.task4_rel1rt(find(d2.group == 1)); d2.task4_rel2rt(find(d2.group == 1))];
% m1 = nanmean(vec); n1 = length(vec)/3; sem1 = 1.96*(nanstd(vec)/sqrt(n1));
% 
% vec = [d2.task4_rel0rt(find(d2.group == 2)); d2.task4_rel1rt(find(d2.group == 2)); d2.task4_rel2rt(find(d2.group == 2))];
% m2 = nanmean(vec); n2 = length(vec)/3; sem2 = 1.96*(nanstd(vec)/sqrt(n2));
% 
% vec = [d2.task4_rel0rt(find(d2.group == 3)); d2.task4_rel1rt(find(d2.group == 3)); d2.task4_rel2rt(find(d2.group == 3))];
% m3 = nanmean(vec); n3 = length(vec)/3; sem3 = 1.96*(nanstd(vec)/sqrt(n3));
% 
% figure(11); hold on;
% b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([2], [m2]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% b = bar([3], [m3]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
% 
% p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
% p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
% p = plot([3 3], [m3+sem3 m3-sem3]); p(1).Color = 'k';
% 
% % plot([0 12], [0.33 0.33], 'k:');
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [0 4];
% xax.TickValues = [1 2 3];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n2)], ['adult, N=' num2str(n3)]};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 6000;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [0 3000 6000];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str((ylim_lo+ylim_hi)/2, '%1.0f'), num2str(ylim_hi, '%1.0f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
% 
% ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
%    
% legend off;
% % legend({'Rel-0', 'Rel-1', 'Rel-2'})
% box off; 
% % legend('box', 'off');
% % legend('location', 'eastoutside');
% pbaspect([2 1 2])
% 
% print(fullfile(rootdir, 'plots', 'devti_relationalreasoning_rt_maineffectofgroup'), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalreasoning_rt_maineffectofgroup'), '-depsc')
% 
% hold off;

%% Task 5: Relational Inference
% Remove data identified as outliers by Schlichting et al., 2017. 
temp = load(fullfile(rootdir, 'taskperf', 'task5_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp;

% % Remove data from subjects with accuracy < .05.
if strcmp(accuracycutoff, 'yes')
    d2.task5_GenDIRrt(d2.task5_GenDIRacc < .5) = NaN;
    d2.task5_GenINDrt(d2.task5_GenINDacc < .5) = NaN;
end

m1 = nanmean(d2.task5_GenDIRrt(find(d2.group == 1))); n1 = length(d2.task5_GenDIRrt(find(d2.group == 1))); sem1 = 1.96*(nanstd(d2.task5_GenDIRrt(find(d2.group == 1)))/sqrt(n1));
m2 = nanmean(d2.task5_GenINDrt(find(d2.group == 1))); n2 = length(d2.task5_GenINDrt(find(d2.group == 1))); sem2 = 1.96*(nanstd(d2.task5_GenINDrt(find(d2.group == 1)))/sqrt(n1));

m3 = nanmean(d2.task5_GenDIRrt(find(d2.group == 2))); n3 = length(d2.task5_GenDIRrt(find(d2.group == 2))); sem3 = 1.96*(nanstd(d2.task5_GenDIRrt(find(d2.group == 2)))/sqrt(n2));
m4 = nanmean(d2.task5_GenINDrt(find(d2.group == 2))); n4 = length(d2.task5_GenINDrt(find(d2.group == 2))); sem4 = 1.96*(nanstd(d2.task5_GenINDrt(find(d2.group == 2)))/sqrt(n2));

m5 = nanmean(d2.task5_GenDIRrt(find(d2.group == 3))); n5 = length(d2.task5_GenDIRrt(find(d2.group == 3))); sem5 = 1.96*(nanstd(d2.task5_GenDIRrt(find(d2.group == 3)))/sqrt(n3));
m6 = nanmean(d2.task5_GenINDrt(find(d2.group == 3))); n6 = length(d2.task5_GenINDrt(find(d2.group == 3))); sem6 = 1.96*(nanstd(d2.task5_GenINDrt(find(d2.group == 3)))/sqrt(n3));

figure(12); hold on;
b = bar([1], [m1]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; barwidth = 0.8; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([2], [m2]); color = [50 180 100]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([4], [m3]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([5], [m4]); color = [50 100 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;
b = bar([7], [m5]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 1;
b = bar([8], [m6]); color = [100 50 180]/255; b.FaceColor = color; b.EdgeColor = 'none'; b.BarWidth = barwidth; b.FaceAlpha = 0.15;

p = plot([1 1], [m1+sem1 m1-sem1]); p(1).Color = 'k';
p = plot([2 2], [m2+sem2 m2-sem2]); p(1).Color = 'k';
p = plot([4 4], [m3+sem3 m3-sem3]); p(1).Color = 'k';
p = plot([5 5], [m4+sem4 m4-sem4]); p(1).Color = 'k';
p = plot([7 7], [m5+sem5 m5-sem5]); p(1).Color = 'k';
p = plot([8 8], [m6+sem6 m6-sem6]); p(1).Color = 'k';

% plot([0 9], [0.5 0.5], 'k:');

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 9];
xax.TickValues = [1.5 4.5 7.5];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {['child, N=' num2str(n1)], ['adol., N=' num2str(n2)], ['adult, N=' num2str(n3)]};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 6000;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [0 3000 6000];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), num2str(3000, '%1.0f'), num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;

ylabel('Reaction Time (ms)', 'FontAngle', 'italic')
   
% legend off;
legend({'Direct', 'Inference'})
box off; 
legend('box', 'off');
legend('location', 'eastoutside');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', 'devti_relationalinference_general_rt'), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', 'devti_relationalinference_general_rt'), '-depsc')

hold off;

[H, P, CI, stats] = ttest(d2.task5_GenDIRrt(find(d2.group == 1)), d2.task5_GenINDrt(find(d2.group == 1)))
[H, P, CI, stats] = ttest(d2.task5_GenDIRrt(find(d2.group == 2)), d2.task5_GenINDrt(find(d2.group == 2)))
[H, P, CI, stats] = ttest(d2.task5_GenDIRrt(find(d2.group == 3)), d2.task5_GenINDrt(find(d2.group == 3)))

[H, P, CI, stats] = ttest2(d2.task5_GenDIRrt(find(d2.group == 1)), d2.task5_GenDIRrt(find(d2.group == 2)))
[H, P, CI, stats] = ttest2(d2.task5_GenDIRrt(find(d2.group == 2)), d2.task5_GenDIRrt(find(d2.group == 3)))

[H, P, CI, stats] = ttest2(d2.task5_GenINDrt(find(d2.group == 2)), d2.task5_GenINDrt(find(d2.group == 3)))

%% Evaluate speed-accuracy trade-off

% Task 1
temp = load(fullfile(rootdir, 'taskperf', 'task1_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp removesubid;

figure(13); 
modelspec = 'assoc ~ assocrt'; 
mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
hold on;
p = plot(mdl);

color = [128 128 128]/255; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = [1 1 1];
p(1).MarkerSize = markersize;
p(2).Color = color;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

legend off;

a = gca;
a.YLabel.String = {'Accuracy (proportion correct)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Reaction Time (ms)'};
a.XLabel.FontSize = fontsize;
% title('AB Direct')
title('');
a.TitleFontSizeMultiplier = 2;

% xaxis
xlimlo = 1000; xlimhi = round(max(d2.task1_rt_dirrt4), -2);
xax = get(gca, 'xaxis');
xax.Limits = [xlimlo xlimhi];
dif = xlimhi - xlimlo;
xax.TickValues = round([xlimlo xlimlo+(dif/4) xlimlo+2*(dif/4) xlimlo+3*(dif/4) xlimhi], -1);
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
legend off;
% legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
box off; 
% legend('box', 'off');
% legend('location', 'southeast');
pbaspect([1 1 1])

[r, p] = corrcoef(d2.task1_acc_dirperf4, d2.task1_rt_dirrt4);
text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.4f')], ...
    'FontSize', fontsize)

print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')

hold off;

figure(14); 
modelspec = 'infer ~ inferrt'; 
mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
hold on;
p = plot(mdl); 

color = [128 128 128]/255; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = [1 1 1];
p(1).MarkerSize = markersize;
p(2).Color = color;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

legend off;

a = gca;
a.YLabel.String = {'Accuracy (proportion correct)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Reaction Time (ms)'};
a.XLabel.FontSize = fontsize;
% title('AC Inference')
title('');
a.TitleFontSizeMultiplier = 2;

% xaxis
xlimlo = 1000; xlimhi = round(max(d2.task1_rt_dirrt4), -2);
xax = get(gca, 'xaxis');
xax.Limits = [xlimlo xlimhi];
dif = xlimhi - xlimlo;
xax.TickValues = round([xlimlo xlimlo+(dif/4) xlimlo+2*(dif/4) xlimlo+3*(dif/4) xlimhi], -1);
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
legend off;
% legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
box off; 
% legend('box', 'off');
% legend('location', 'southeast');
pbaspect([1 1 1])

[r, p] = corrcoef(d2.task1_acc_ACperf, d2.task1_rt_ACrt)
text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.3f')], ...
    'FontSize', fontsize)

print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')

hold off;

% % Task 4
% temp = load(fullfile(rootdir, 'taskperf', 'task4_behavior_n=78.mat'));
% removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
% removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
% removesubid = unique([removetemp1 removetemp2]);
% clear temp removetemp1 removetemp2;
% d2 = d(~ismember(d.subID, removesubid), :); clear temp removesubid;
% 
% figure(15); 
% modelspec = 'task4_rel0acc ~ task4_rel0rt'; 
% mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
% hold on;
% p = plot(mdl); 
% 
% color = [128 128 128]/255; 
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = [1 1 1];
% p(1).MarkerSize = markersize;
% p(2).Color = color;
% p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;
% 
% legend off;
% 
% a = gca;
% a.YLabel.String = {'Accuracy (proportion correct)'};
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% a.XLabel.String = {'Reaction Time (ms)'};
% a.XLabel.FontSize = fontsize;
% % title('Subregion: CA1')
% title('Rel-0')
% a.TitleFontSizeMultiplier = 2;
% 
% % xaxis
% xlimlo = 1000; xlimhi = 8000;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% dif = xlimhi - xlimlo;
% xax.TickValues = round([xlimlo xlimlo+(dif/4) xlimlo+2*(dif/4) xlimlo+3*(dif/4) xlimhi], -1);
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 1;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
%    
% legend off;
% % legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
% box off; 
% % legend('box', 'off');
% % legend('location', 'southeast');
% pbaspect([1 1 1])
% 
% [r, p] = corrcoef(d2.task4_rel0acc, d2.task4_rel0rt)
% text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.3f')], ...
%     'FontSize', fontsize)
% 
% print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')
% 
% hold off;
% 
% figure(16); 
% modelspec = 'task4_rel1acc ~ task4_rel1rt'; 
% mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
% hold on;
% p = plot(mdl); 
% 
% color = [128 128 128]/255; 
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = [1 1 1];
% p(1).MarkerSize = markersize;
% p(2).Color = color;
% p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;
% 
% legend off;
% 
% a = gca;
% a.YLabel.String = {'Accuracy (proportion correct)'};
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% a.XLabel.String = {'Reaction Time (ms)'};
% a.XLabel.FontSize = fontsize;
% % title('Subregion: CA1')
% title('Rel-1')
% a.TitleFontSizeMultiplier = 2;
% 
% % xaxis
% xlimlo = 1000; xlimhi = 8000;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% dif = xlimhi - xlimlo;
% xax.TickValues = round([xlimlo xlimlo+(dif/4) xlimlo+2*(dif/4) xlimlo+3*(dif/4) xlimhi], -1);
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 1;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
%    
% legend off;
% % legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
% box off; 
% % legend('box', 'off');
% % legend('location', 'southeast');
% pbaspect([1 1 1])
% 
% [r, p] = corrcoef(d2.task4_rel1acc, d2.task4_rel1rt)
% text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.3f')], ...
%     'FontSize', fontsize)
% 
% print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')
% 
% hold off;
% 
% figure(17); 
% modelspec = 'task4_rel2acc ~ task4_rel2rt'; 
% mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
% hold on;
% p = plot(mdl); 
% 
% color = [128 128 128]/255; 
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = [1 1 1];
% p(1).MarkerSize = markersize;
% p(2).Color = color;
% p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;
% 
% legend off;
% 
% a = gca;
% a.YLabel.String = {'Accuracy (proportion correct)'};
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% a.XLabel.String = {'Reaction Time (ms)'};
% a.XLabel.FontSize = fontsize;
% % title('Subregion: CA1')
% title('Rel-2')
% a.TitleFontSizeMultiplier = 2;
% 
% % xaxis
% xlimlo = 1000; xlimhi = 8000;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% dif = xlimhi - xlimlo;
% xax.TickValues = round([xlimlo xlimlo+(dif/4) xlimlo+2*(dif/4) xlimlo+3*(dif/4) xlimhi], -1);
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
% ylim_lo = 0; ylim_hi = 1;
% yax = get(gca,'yaxis');
% yax.Limits = [ylim_lo ylim_hi];
% yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
% yax.FontAngle = fontangle;
%    
% legend off;
% % legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
% box off; 
% % legend('box', 'off');
% % legend('location', 'southeast');
% pbaspect([1 1 1])
% 
% [r, p] = corrcoef(d2.task4_rel2acc, d2.task4_rel2rt)
% text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.3f')], ...
%     'FontSize', fontsize)
% 
% print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')
% 
% hold off;

% Task 5
temp = load(fullfile(rootdir, 'taskperf', 'task5_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp removesubid;

figure(18); 
modelspec = 'task5_GenDIRacc ~ task5_GenDIRrt'; 
mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
hold on;
p = plot(mdl); 

color = [128 128 128]/255; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = [1 1 1];
p(1).MarkerSize = markersize;
p(2).Color = color;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

legend off;

a = gca;
a.YLabel.String = {'Accuracy (proportion correct)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Reaction Time (ms)'};
a.XLabel.FontSize = fontsize;
% title('Gen Direct')
title('');
a.TitleFontSizeMultiplier = 2;

% xaxis
xlimlo = 0; xlimhi = 5000;
xax = get(gca, 'xaxis');
xax.Limits = [xlimlo xlimhi];
xax.TickValues = [xlimlo xlimhi/4 2*(xlimhi/4) 3*(xlimhi/4) xlimhi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
legend off;
% legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
box off; 
% legend('box', 'off');
% legend('location', 'southeast');
pbaspect([1 1 1])

[r, p] = corrcoef(d2.task5_GenDIRacc, d2.task5_GenDIRrt)
text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.3f')], ...
    'FontSize', fontsize)

print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')

hold off;

figure(19); 
modelspec = 'task5_GenINDacc ~ task5_GenINDrt'; 
mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
hold on;
p = plot(mdl); 

color = [128 128 128]/255; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = [1 1 1];
p(1).MarkerSize = markersize;
p(2).Color = color;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

legend off;

a = gca;
a.YLabel.String = {'Accuracy (proportion correct)'};
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Reaction Time (ms)'};
a.XLabel.FontSize = fontsize;
% title('Gen Inference')
title('');
a.TitleFontSizeMultiplier = 2;

% xaxis
xlimlo = 0; xlimhi = round(max(d2.task1_rt_dirrt4), -2);
xax = get(gca, 'xaxis');
xax.Limits = [xlimlo xlimhi];
xax.TickValues = [xlimlo xlimhi/4 2*(xlimhi/4) 3*(xlimhi/4) xlimhi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {num2str(xlimlo, '%1.0f'), num2str(xlimhi/4, '%1.0f'), num2str(2*(xlimhi/4), '%1.0f'), num2str(3*(xlimhi/4),'%1.0f'), num2str(xlimhi, '%1.0f')};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
ylim_lo = 0; ylim_hi = 1;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
legend off;
% legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
box off; 
% legend('box', 'off');
% legend('location', 'southeast');
pbaspect([1 1 1])

[r, p] = corrcoef(d2.task5_GenINDacc, d2.task5_GenINDrt)
text(xlimlo+100, ylim_lo+0.10, ['r = ' num2str(r(1, 2), '%1.3f') ', p = ' num2str(p(1, 2), '%1.3f')], ...
    'FontSize', fontsize)

print(fullfile(rootdir, 'plots', ['speedaccuracytradeoff_' modelspec]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['speedaccuracytradeoff_' modelspec]), '-depsc')

hold off;

%% Export data for ANOVA in SPSS

% All data.
filename = sprintf('devti_data_beh_forSPSS_%s', datestr(now,'yyyymmdd'));
save(fullfile(rootdir, 'supportFiles', filename), 'd')
writetable(d, fullfile(rootdir, 'supportFiles', [filename '.csv']))

%% Task 1 and Task 3: Final Direct and Inference, Schlichting et al., 2017
% Remove data identified as outliers by Schlichting et al., 2017. 
temp = load(fullfile(rootdir, 'taskperf', 'task1_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp removesubid;
filename = sprintf('devti_data_beh_forSPSS_task1outliersremoved_%s', datestr(now,'yyyymmdd'));
save(fullfile(rootdir, 'supportFiles', filename), 'd2')
writetable(d2, fullfile(rootdir, 'supportFiles', [filename '.csv']))
clear d2;

% %% Task 4: Relational reasoning, Schlichting et al., 2017, Crone and Bunge
% % Remove data identified as outliers by Schlichting et al., 2017. 
% temp = load(fullfile(rootdir, 'taskperf', 'task4_behavior_n=78.mat'));
% removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
% removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
% removesubid = unique([removetemp1 removetemp2]);
% clear temp removetemp1 removetemp2;
% d2 = d(~ismember(d.subID, removesubid), :); clear temp;
% filename = sprintf('devti_data_beh_forSPSS_task4outliersremoved_%s', datestr(now,'yyyymmdd'));
% save(fullfile(rootdir, 'supportFiles', filename), 'd2')
% writetable(d2, fullfile(rootdir, 'supportFiles', [filename '.csv']))
% clear d2;

%% Task 5: Relational Inference, Wendelken and Bunge
% Remove data identified as outliers by Schlichting et al., 2017. 
temp = load(fullfile(rootdir, 'taskperf', 'task5_behavior_n=78.mat'));
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
removesubid = unique([removetemp1 removetemp2]);
clear temp removetemp1 removetemp2;
d2 = d(~ismember(d.subID, removesubid), :); clear temp;
filename = sprintf('devti_data_beh_forSPSS_task5outliersremoved_%s', datestr(now,'yyyymmdd'));
save(fullfile(rootdir, 'supportFiles', filename), 'd2')
writetable(d2, fullfile(rootdir, 'supportFiles', [filename '.csv']))
clear d2;
