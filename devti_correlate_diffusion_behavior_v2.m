% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

% Note that includeBoolPerf and includeBoolRT 
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

% Adjust what measure is used to quantify assoc and infer. Default is task1_acc_dirperf4 and task1_acc_ACperf.
d.assoc = mean([d.task1_acc_dirperf1 d.task1_acc_dirperf2 d.task1_acc_dirperf3 d.task1_acc_dirperf4], 2);

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

yc_color  = [220 220 220]/255; %[50 180 100]/255; 
oc_color = [169 169 169]/255; %[50 100 180]/255; 
a_color = [105 105 105]/255; %[100 50 180]/255;

ifof = [142 198 255]/255; % light blue
tpc =  [27 102 87]/255; % dark turquoise
mdlfspl = [42, 102, 0]/255; % green
mdlfang =  [72 146 73]/255; % medium sea green

% Rename variables to be consistent with manuscript.
d = renamevars(d, ["assoc", "task5_GenDIRacc", "infer", "task5_GenINDacc"], ["Mmatching", "Pmatching", "Minference", "Pinference"]);

%% Models.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear main effects model: CA1.
modelspec = 'ca1 ~ Minference + Pinference + age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: CA1.
modelspec = 'ca1 ~ Minference*age + Pinference*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: CA1.
modelspec = 'ca1 ~ Mmatching + Pmatching + age^2 + iq';  
% modelspec = 'ca1 ~ Mmatching*age + Pmatching*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: CA1.
modelspec = 'ca1 ~ Mmatching*age + Pmatching*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear main effects model: CA23.
modelspec = 'ca23 ~ Minference + Pinference + sex*age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: CA23.
modelspec = 'ca23 ~ Minference*age + Pinference*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: CA23.
modelspec = 'ca23 ~ Mmatching + Pmatching + sex*age^2 + iq';  
%modelspec = 'ca23 ~ Mmatching*age + Pmatching*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: CA23.
modelspec = 'ca23 ~ Mmatching*age + Pmatching*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear main effects model: DG.
modelspec = 'dg ~ Minference + Pinference + age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: DG.
modelspec = 'dg ~ Minference*age + Pinference*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: DG.
modelspec = 'dg ~ Mmatching + Pmatching + age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: DG.
modelspec = 'dg ~ Mmatching*age + Pmatching*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear main effects model: SUB.
modelspec = 'sub ~ Minference + Pinference + age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: SUB.
modelspec = 'sub ~ Minference*age + Pinference*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: SUB.
modelspec = 'sub ~ Mmatching + Pmatching + age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Linear main effects model: SUB.
modelspec = 'sub ~ Mmatching*age + Pmatching*age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT: Gen Direct predicts CA1 microstructure.

% Linear main effects model.
modelspec = 'ca1 ~ Mmatching + Pmatching + age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

% Plot.
figure(1);
hold on;
p = plotAdjustedResponse(mdl, 'Pmatching'); color = ca1; 
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

a.XLabel.String = {'Gen Direct Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
title('Cornu Ammonis 1 (CA1)')

if strcmp(wm, 'fa')
    ylim_lo = 0.10; ylim_hi = 0.30;
elseif strcmp(wm, 'md')
    ylim_lo = 0; ylim_hi = 1.2;  
end

% xaxis
xlim_lo = 0.88; xlim_hi = 1.00;
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

%% PLOT: AC Inference predicts DG microstructure.

% Linear main effects model.
modelspec = 'dg ~ Minference + Pinference + age^2 + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% % Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

close all;
f = figure('visible', 'on');
set(gcf,'Visible','on');              
set(0,'DefaultFigureVisible','on');

figure(2); hold on;
p = plotAdjustedResponse(mdl, 'Minference'); color = ca1; 
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

%% PLOT: GEN Inference predicts SUB microstructure.

% Linear main effects model.
modelspec = 'sub ~ Minference + Pinference + age + iq';  

% Remove outliers.
removeidx = findoutliers(d, modelspec);
d2  = d; d2(removeidx, :) = []; 

% % Fit regression model.
mdl = fitlm(d2, modelspec)
n = size(d2, 1);
disp(['N = ' num2str(n)])
disp(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
disp('');

close all;
f = figure('visible', 'on');
set(gcf,'Visible','on');              
set(0,'DefaultFigureVisible','on');

figure(3); hold on;
p = plotAdjustedResponse(mdl, 'Pinference'); color = ca1; 
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

a.XLabel.String = {'Gen Inference Accuracy'; '(adjusted)'};
a.XLabel.FontSize = fontsize;
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

n = size(d2, 1);
print(fullfile(rootdir, 'plots', ['adjusted_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['adjusted_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off; 

% Get correlations among variables of interest. NOTE: d2 variables have all
% been demeaned above, except age.
voi = [d2.age d2.head d2.body d2.tail d2.ca1 d2.ca23 d2.dg d2.sub d2.Mmatching d2.Minference d2.Pmatching d2.Pinference];
[rval, pval] = corr(voi);




