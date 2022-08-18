% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

accuracycutoff = 'no'; % selecting "no" going forward (7/8/2022) because includeBoolPerf and includeBookRT 
% exclude subjects who had poor performance due to not understanding the task, they had a complex way of 
% determining whether poor performance was a lack of understanding the task or just poor performance

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Load behavioral data.
load(fullfile(rootdir, 'supportFiles', 'devti_data_beh_20220705.mat'));
iq = readtable(fullfile(rootdir, 'supportFiles', 'devti_iq_all.xlsx')); 
iq.Properties.VariableNames(1) = {'subID'}; iq.Properties.VariableNames(2) = {'iq'};

% Load diffusion data.
wm = 'fa';
load(fullfile(rootdir, 'supportFiles', ['devti_data_' wm '_20220705.mat']));
mri = m(:, [1 4:end]); clear m;
% Scale md values for analysis and visualization.
if strcmp(wm, 'md')
    mri(:, 5:end) = array2table(table2array(mri(:, 5:end)).*1000);
end

% Rectify subIDs.
idx = find(ismember(mri.subID, beh.subID));
mri = mri(idx, :); clear idx;
idx = find(ismember(beh.subID, mri.subID));
beh = beh(idx, :); clear idx;
idx = find(ismember(iq.subID, beh.subID));
iq = iq(idx, :); clear idx;

% Concatenate beh and mri data into one table.
d = [beh(:, 1:3) mri(:, 2:4) beh(:, 4:end) iq(:, 2) mri(:, 5:end)];

% Remove the random spaces and other issues in the column names.
d.Properties.VariableNames = strrep(d.Properties.VariableNames, ' ', '');
d.Properties.VariableNames = strrep(d.Properties.VariableNames, '|', '');
d.Properties.VariableNames = strrep(d.Properties.VariableNames, '&', '');

% Adjust what measure is used to quantify assoc and infer. Default is
% task1_acc_dirperf4 and task1_acc_ACperf.
% d.assoc = mean([d.task1_acc_dirperf1 d.task1_acc_dirperf2 d.task1_acc_dirperf3 d.task1_acc_dirperf4], 2);
% d.infer = d.task1_acc_ACperfABCBCC;

% Tell Matlab that sex and age group are categorical variables.
d.sex = categorical(d.sex);
% d.group = categorical(d.group);
clear beh mri;

% Add anterior and posterior.
d.anterior = d.head;
d.posterior = mean([d.body d.tail], 2);

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

ifof = [142 198 255]/255; % light blue
tpc =  [27 102 87]/255; % dark turquoise
mdlfspl = [42, 102, 0]/255; % green
mdlfang =  [72 146 73]/255; % medium sea green

%% Inference with memory, AB Direct and AC Inference tasks, Schlichting et al., 2017

% Linear main effects model.
modelspec = 'dg ~ assoc + infer + task5_GenDIRacc + task5_GenINDacc + age^2 + iq';  
% modelspec = 'infer ~ ca1 + ca23 + dg + sub + age + assoc';  

% Rmeoving multivariate outliers for 'task5_GenDIRacc ~ ca1 + ca23 + dg + sub + age + sex + iq'leads to all 1s
% for Gen Diracc.

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% % Identify subjects that were below chance on either task.
% idx1 = d2.subID(find(d2.assoc <= . | d2.infer <= .33));
% idx4 = d2.subID(find(d2.task4_rel0acc <= .33 | d2.task4_rel1acc <= .33 | d2.task4_rel2acc <= .33));
% idx5 = d2.subID(find(d2.task5_GenDIRacc <= .5 | d2.task5_GenINDacc <= .50));

if strcmp(accuracycutoff, 'yes')
    removesub = idx1; % idx1 for inference memory (object task), idx5 for inference no memory (ball task)
    d2 = d2(~ismember(d2.subID, removesub), :);
end

% % Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');
% p = plotAdjustedResponse(mdl, 'ca1'); color = ca1; 
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = color;
% p(1).MarkerSize = markersize;
% p(2).Color = color;
% p(2).LineWidth = linewidth;
% % p(3).LineStyle = 'none';
% % p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; %CI = (p(4).YData - p(3).YData)/2; 
% %patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;
figure(1);
hold on;
p = plotAdjustedResponse(mdl, 'infer'); color = ca1; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = [128 128 128]/255;
p(2).LineWidth = linewidth;
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, d2.group, [yc_color; oc_color; a_color]);
marker = 'o';
s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;
% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], [128 128 128]/255, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none';
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD) (adjusted)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'AC Inference Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA1')
title('Dentate Gyrus (DG)')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xlim_lo = 0; xlim_hi = 1.00;
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.TickLabels = {num2str(xlim_lo, '%1.2f'), num2str((xlim_lo+xlim_hi)/2, '%1.2f'), num2str(xlim_hi, '%1.2f')};
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
% legend off;
legend({'', 'fit', 'child', 'adol', 'adult'})
box off; 
legend('box', 'off');
legend('location', 'northeast');
pbaspect([1 1 1])

n = size(d2, 1);
print(fullfile(rootdir, 'plots', ['adjusted_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['adjusted_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off; 

figure(2);
hold on;
p = plotAdjustedResponse(mdl, 'assoc'); color = ca1; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = [128 128 128]/255;
p(2).LineWidth = linewidth;
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, d2.group, [yc_color; oc_color; a_color]);
marker = 'o';
s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;
% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], [128 128 128]/255, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none';
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD) (adjusted)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'AB Direct Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA1')
title('Dentate Gyrus (DG)')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xlim_lo = 0; xlim_hi = 1.00;
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.TickLabels = {num2str(xlim_lo, '%1.2f'), num2str((xlim_lo+xlim_hi)/2, '%1.2f'), num2str(xlim_hi, '%1.2f')};
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
% legend off;
legend({'', 'fit', 'child', 'adol', 'adult'})
box off; 
legend('box', 'off');
legend('location', 'northeast');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['adjusted_' wm '_' modelspec '_n=' num2str(n) '_2']), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['adjusted_' wm '_' modelspec '_n=' num2str(n) '_2']), '-depsc')

% Mean center continuous variables for modelling.
d2(:, 5:end) = array2table(double(table2array(d2(:, 5:end)) - nanmean(table2array(d2(:, 5:end)), 1)));
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

hold off; clear d2;

%% Inference without memory, ball task

modelspec = 'sub ~ assoc + infer + task5_GenDIRacc + task5_GenINDacc + age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% % Identify subjects that were below chance on either task.
% idx5 = d.subID(find(d.task5_GenDIRacc <= .5 | d.task5_GenINDacc <= .50));
% removesub = idx5; % idx1 for inference memory (object task), idx5 for inference no memory (ball task)
% d2 = d(~ismember(d.subID, removesub), :);

% Mean center continuous variables for modelling.
% d2(:, 5:end) = array2table(double(table2array(d2(:, 5:end)) - nanmean(table2array(d2(:, 5:end)), 1)));

figure(3);
% Fit regression model.
mdl = fitlm(d2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
% p = plotAdjustedResponse(mdl, 'ca1'); color = ca1; 
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = color;
% p(1).MarkerSize = markersize;
% p(2).Color = color;
% p(2).LineWidth = linewidth;
% % p(3).LineStyle = 'none';
% % p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; %CI = (p(4).YData - p(3).YData)/2; 
% %patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;
figure(4); hold on;
p = plotAdjustedResponse(mdl, 'task5_GenINDacc'); color = dg; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = [128 128 128]/255;
p(2).LineWidth = linewidth;
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, d2.group, [yc_color; oc_color; a_color]);
marker = 'o';
s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;
% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], [128 128 128]/255, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none';
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.XLabel.String = {'Mean Diffusivity (MD) (adjusted)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Gen Inference Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA1')
title('Subiculum (SUB)')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xlim_lo = 0; xlim_hi = 1.00;
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.TickLabels = {num2str(xlim_lo, '%1.2f'), num2str((xlim_lo+xlim_hi)/2, '%1.2f'), num2str(xlim_hi, '%1.2f')};
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
% legend off;
legend({'', 'fit', 'child', 'adol', 'adult'})
box off; 
legend('box', 'off');
legend('location', 'southeast');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['adjusted_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['adjusted_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

figure(5); figure(6); hold off
figure(8); hold on;
p = plotAdjustedResponse(mdl, 'task5_GenDIRacc'); color = dg; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = [128 128 128]/255;
p(2).LineWidth = linewidth;
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, d2.group, [yc_color; oc_color; a_color]);
marker = 'o';
s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;
% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], [128 128 128]/255, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none';
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

% Mean center continuous variables for modelling.
d2(:, 5:end) = array2table(double(table2array(d2(:, 5:end)) - nanmean(table2array(d2(:, 5:end)), 1)));
mdl = fitlm(d2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
n = size(d2, 1);
% disp('CA1')
display(['adjr2 = ' num2str(mdl.Rsquared.Adjusted)])
display(['rmse = ' num2str(mdl.RMSE)])
% display(['beta = ' num2str(mdl.Coefficients.Estimate(2))]); %NOTE: This will change with the number of predictors.
% display(['p (for beta) = ' num2str(mdl.Coefficients.pValue(2))]); %NOTE: This will change with the number of predictors.
display(['n = ' num2str(n)]); 
% clear d2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.XLabel.String = {'Mean Diffusivity (MD) (adjusted)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Gen Direct Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA1')
title('Subiculum (SUB)')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xlim_lo = 0; xlim_hi = 1.00;
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xax.TickLabels = {num2str(xlim_lo, '%1.2f'), num2str((xlim_lo+xlim_hi)/2, '%1.2f'), num2str(xlim_hi, '%1.2f')};
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
% legend off;
legend({'', 'fit', 'child', 'adol', 'adult'})
box off; 
legend('box', 'off');
legend('location', 'northwest');
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['adjusted_' wm '_' modelspec '_n=' num2str(n) '_2']), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['adjusted_' wm '_' modelspec '_n=' num2str(n) '_2']), '-depsc')

hold off;

% Mean center continuous variables for modelling.
d2(:, 5:end) = array2table(double(table2array(d2(:, 5:end)) - nanmean(table2array(d2(:, 5:end)), 1)));
mdl = fitlm(d2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
n = size(d2, 1);
% disp('CA1')
display(['adjr2 = ' num2str(mdl.Rsquared.Adjusted)])
display(['rmse = ' num2str(mdl.RMSE)])
% display(['beta = ' num2str(mdl.Coefficients.Estimate(2))]); %NOTE: This will change with the number of predictors.
% display(['p (for beta) = ' num2str(mdl.Coefficients.pValue(2))]); %NOTE: This will change with the number of predictors.
display(['n = ' num2str(n)]); 
% clear d2;

%% Age AC Inference interaction for subiculum: 
% Inference with memory, AB Direct and AC Inference tasks, Schlichting et al., 2017

% Linear main effects model.
modelspec = 'sub ~ assoc*age + infer*age + task5_GenDIRacc*age + task5_GenINDacc*age + iq';  
% modelspec = 'infer ~ ca1 + ca23 + dg + sub + age + assoc';  

% Rmeoving multivariate outliers for 'task5_GenDIRacc ~ ca1 + ca23 + dg + sub + age + sex + iq'leads to all 1s
% for Gen Diracc.

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% % Identify subjects that were below chance on either task.
% idx1 = d2.subID(find(d2.assoc <= . | d2.infer <= .33));
% idx4 = d2.subID(find(d2.task4_rel0acc <= .33 | d2.task4_rel1acc <= .33 | d2.task4_rel2acc <= .33));
% idx5 = d2.subID(find(d2.task5_GenDIRacc <= .5 | d2.task5_GenINDacc <= .50));

if strcmp(accuracycutoff, 'yes')
    removesub = idx1; % idx1 for inference memory (object task), idx5 for inference no memory (ball task)
    d2 = d2(~ismember(d2.subID, removesub), :);
end

% % Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');
% p = plotAdjustedResponse(mdl, 'ca1'); color = ca1; 
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = color;
% p(1).MarkerSize = markersize;
% p(2).Color = color;
% p(2).LineWidth = linewidth;
% % p(3).LineStyle = 'none';
% % p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; %CI = (p(4).YData - p(3).YData)/2; 
% %patch([x fliplr(x)], [y-CI fliplr(y+CI)], color, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;
mdl = fitlm(d2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
p = plotAdjustedResponse(mdl, 'age');
hold on;

color = [128 128 128]/255;
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = [128 128 128]/255;
p(2).LineWidth = linewidth;
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, d2.group, [yc_color; oc_color; a_color]);
marker = 'o';
s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;

% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], [128 128 128]/255, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none';
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

% % Mean center continuous variables for modelling.
% m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
% mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
% n = mdl.NumObservations;
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
% clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'Age'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA23')
title('')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [5 30];
xax.TickValues = [5 10 15 20 25 30];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'5', '10', '15', '20', '25', '30'};
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
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
legend off;
% legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
box off; 
% legend('box', 'off');
% legend('location', 'southeast');
pbaspect([1 1 1]);
title('');

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n) '_age']), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n) '_age']), '-depsc')

hold off;

figure(13);
p = plotAdjustedResponse(mdl, 'task5_GenINDacc');
hold on;

color = [128 128 128]/255;
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = [128 128 128]/255;
p(2).LineWidth = linewidth;
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, d2.group, [yc_color; oc_color; a_color]);
marker = 'o';
s(1).MarkerSize = markersize; s(1).Marker = marker; s(1).MarkerEdgeColor = 'w'; s(1).MarkerFaceColor = yc_color;
s(2).MarkerSize = markersize; s(2).Marker = marker; s(2).MarkerEdgeColor = 'w'; s(2).MarkerFaceColor = oc_color;
s(3).MarkerSize = markersize; s(3).Marker = marker; s(3).MarkerEdgeColor = 'w'; s(3).MarkerFaceColor = a_color;

% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], [128 128 128]/255, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none';
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

% % Mean center continuous variables for modelling.
% m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
% mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
% n = mdl.NumObservations;
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
% clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'Age'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA23')
title('')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 1];
xax.TickValues = [0 0.5 1.0];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'0.00', '0.50' '1.00'};
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
yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.2f'), num2str(ylim_hi, '%1.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
legend off;
% legend({'DG', '', 'CA1', '', 'CA23', '', 'SUB', ''})
box off; 
% legend('box', 'off');
% legend('location', 'southeast');
pbaspect([1 1 1]);
title('');

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n) '_task5_GenINDacc']), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n) '_task5_GenINDacc']), '-depsc')

hold off;

figure(14)
p = plotInteraction(mdl, 'age', 'task5_GenINDacc', 'predictions');
p(1).Color = yc_color; p(1).LineStyle = '-'; p(1).LineWidth = 8;
p(2).Color = oc_color; p(2).LineStyle = '-'; p(2).LineWidth = 8;
p(3).Color = a_color; p(3).LineStyle = '-'; p(3).LineWidth = 8;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'Gen Inference Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA23')
title('')

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 1];
xax.TickValues = [0 .5 1];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'0.00', '0.05', '1.00'};
xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0.10 0.30];
yax.TickValues = [0.10 0.20 0.30];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {'0.10', '0.20', '0.30'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
% legend off;
legend({'Age', 'Younger (6 year', 'Middle' , 'Older'})
box off; 
legend('box', 'off');
l = legend('location', 'southeast'); l.FontSize = fontsize;
pbaspect([1 1 1]);
title('');

print(fullfile(rootdir, 'plots', ['adjusted_interaction_' wm '_' modelspec '_interaction_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['adjusted_interaction_' wm '_' modelspec '_interaction_n=' num2str(n)]), '-depsc')

hold off;
% Get correlations among variables of interest. NOTE: d2 variables have all
% been demeaned above, except age.
voi = [d2.age d2.head d2.body d2.tail d2.ca1 d2.ca23 d2.dg d2.sub d2.assoc d2.infer d2.task5_GenDIRacc d2.task5_GenINDacc];
[rval, pval] = corr(voi);



