% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

% Note bl-repository.
blprojectid = '5e5672430f7fa65e1d3c9621';

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Select measure and roi.
wm = 'fa'; %'ad', 'rd', 'md'

% Load data.
date = '20211110';
load(fullfile(rootdir, ['supportFiles/devti_data_' wm '_' date '.mat']))

% Tell Matlab that sex and age group are categorical variables.
m.sex = categorical(m.sex);
m.age_group = categorical(m.age_group);

% Scale md values for analysis and visualization.
if strcmp(wm, 'md')
    m(:, 5:end) = array2table(table2array(m(:, 5:end)).*1000);
end
% 
% Set up plot and measure-specific details.
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

gray = [128 128 128]/255;
head = [128 128 128]/255; %[204 204 0]/255; % light burnt yellow
body = [128 128 128]/255; %[204 190 0]/255; % burnt yellow
tail =  [128 128 128]/255; %[210 43 43]/255; % cadmium red

hip = [128 128 128]/255; %[2 129 129]/255; % dark turquoise

ca1color = [128 128 128]/255; %[0 150 255]/255;
ca23color = [128 128 128]/255; %[199 21 133]/255;
dgcolor = [128 128 128]/255; %[75 0 130]/255;
subcolor = [128 128 128]/255; %[60 179 113]/255;

% yc_color  = [50 180 100]/255; 
% oc_color = [50 100 180]/255; 
% a_color = [100 50 180]/255;

if strcmp(wm, 'fa')
    ylim_lo = 0.11; ylim_hi = 0.34;
elseif strcmp(wm, 'md')
    ylim_lo = .7; ylim_hi = 1.4;  
end


%% Hip

% Linear main effects model.
modelspec = 'hip ~ sex';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 

% figure(1);
% Fit regression model.
% mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
% p = plot(mdl);
% hold on;
% 
% color = [0 127 255]/255; % dark blue
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = color;
% p(1).MarkerSize = markersize;
% p(2).Color = gray;
% p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;
% 
% a = gca;
% if strcmp(wm, 'fa')
%     a.YLabel.String = {'Fractional Anisotropy (FA)'};
% elseif strcmp(wm, 'md')
%     a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
% end
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% a.XLabel.String = {'Age'};
% a.XLabel.FontSize = fontsize;
% % title('Subregion: CA1')
% title('')
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [5 30];
% xax.TickValues = [5 10 15 20 25 30];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {'5', '10', '15', '20', '25', '30'};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
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
% print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')
% 
% hold off;

%% Head

% Linear main effects model.
modelspec = 'head ~ sex';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 

figure(1);
% Fit regression model.
mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
p = plot(mdl);
hold on;

color = head; 
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = color;
p(1).MarkerSize = markersize;
p(2).Color = gray;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Sex'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA1')
title('')

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 3];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Female', 'Male'};
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
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

%% Body

% Linear main effects model.
modelspec = 'body ~ sex';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 

% figure(2);
% Fit regression model.
% mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
% p = plot(mdl);
% hold on;
% 
% color = [0 127 255]/255; % dark blue
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = color;
% p(1).MarkerSize = markersize;
% p(2).Color = gray;
% p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;
% 
% a = gca;
% if strcmp(wm, 'fa')
%     a.YLabel.String = {'Fractional Anisotropy (FA)'};
% elseif strcmp(wm, 'md')
%     a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
% end
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% a.XLabel.String = {'Age'};
% a.XLabel.FontSize = fontsize;
% % title('Subregion: CA1')
% title('')
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [5 30];
% xax.TickValues = [5 10 15 20 25 30];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {'5', '10', '15', '20', '25', '30'};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
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
% print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')
% 
% hold off;

%% Tail

% Linear main effects model.
modelspec = 'tail ~ sex';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 

% figure(4);
% % Fit regression model.
% mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% % display(['N = ' num2str(n)])
% % display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
% p = plot(mdl);
% hold on;
% 
% color = [0 127 255]/255; % dark blue
% p(1).Color = color;
% p(1).Marker = 'o';
% p(1).MarkerFaceColor = color;
% p(1).MarkerEdgeColor = color;
% p(1).MarkerSize = markersize;
% p(2).Color = gray;
% p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
% clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;
% 
% a = gca;
% if strcmp(wm, 'fa')
%     a.YLabel.String = {'Fractional Anisotropy (FA)'};
% elseif strcmp(wm, 'md')
%     a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
% end
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% a.XLabel.String = {'Age'};
% a.XLabel.FontSize = fontsize;
% % title('Subregion: CA1')
% title('')
% 
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [5 30];
% xax.TickValues = [5 10 15 20 25 30];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xlabels = {'5', '10', '15', '20', '25', '30'};
% xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% 
% % yaxis
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
% print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
% print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')
% 
% hold off;

%% CA1

% Linear main effects model.
modelspec = 'ca1 ~ age_cov^2';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 

figure(6);
% Fit regression model.
mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
p = plot(mdl);
hold on;

color = ca1color; %[0 127 255]/255; % dark blue
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = color;
p(1).MarkerSize = markersize;
p(2).Color = gray;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Age'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA1')
title('')

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
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

%% CA23

% Linear main effects model.
modelspec = 'ca23 ~ age_cov^2 * sex';   

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 
    
figure(7);
% Fit regression model.
mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
% display(['N = ' num2str(n)])
% display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
p = plotAdjustedResponse(mdl, 'age_cov');
hold on;

color = ca23color; %[240 128 128]/255; % light coral
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = color;
p(1).MarkerSize = markersize;
p(2).Color = gray;
p(2).LineWidth = linewidth;
% p(3).LineStyle = 'none';
% p(4).LineStyle = 'none';
% x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
% patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')

% Get data for plotting the confidence intervals and add CI to plot.
t = array2table(cat(2, p(1).XData', p(1).YData')); t.Properties.VariableNames =  {'x', 'y'};
mdlci = fitlm(t, 'y~x^2');
pci = plot(mdlci);
x = pci(2).XData; y = pci(2).YData; CI = (pci(4).YData - pci(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
pci(1).Marker = 'none'; pci(2).Color = color;
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p pci;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'Age'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA23')
title('')

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

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

figure(8);
p = plotInteraction(mdl, 'sex', 'age_cov', 'predictions');
p(1).Color = color; p(1).LineStyle = '-'; p(1).LineWidth = 10;
p(2).Color = color; p(2).LineStyle = ':'; p(2).LineWidth = 10;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'; '(adjusted)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;

a.XLabel.String = {'Age'};
a.XLabel.FontSize = fontsize;
% title('Subregion: CA23')
title('')

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
yax.Limits = [-0.06 0.10];
yax.TickValues = [-0.06 -0.02 0.02 0.06 0.10];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {'-0.06', '-0.02', '0.02', '0.06', '0.10'};
yax.FontName = fontname;
yax.FontSize = fontsize;
yax.FontAngle = fontangle;
   
% legend off;
legend({'', 'Female', 'Male'})
box off; 
legend('box', 'off');
l = legend('location', 'northwest'); l.FontSize = fontsize;
pbaspect([1 1 1]);
title('');

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_interaction_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_interaction_n=' num2str(n)]), '-depsc')

hold off;

%% DG

% Linear main effects model.
modelspec = 'dg ~ age_cov^2';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 
    
figure(9);
% Fit regression model.
mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
p = plot(mdl);
hold on;

color = dgcolor; %[147 112 219]/255; % medium purple
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = color;
p(1).MarkerSize = markersize;
p(2).Color = gray;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Age'};
a.XLabel.FontSize = fontsize;
% title('Subregion: DG')
title('')

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
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

%% SUB

% Linear main effects model.
modelspec = 'sub ~ age_cov';  

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 
    
figure(10);
% Fit regression model.
mdl = fitlm(m2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
p = plot(mdl);
hold on;

color = subcolor; %[60 179 113]/255; %[64 224 208]/255; % turquoise
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = color;
p(1).MarkerEdgeColor = color;
p(1).MarkerSize = markersize;
p(2).Color = gray;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
x = p(2).XData; y = p(2).YData; CI = (p(4).YData - p(3).YData)/2; 
patch([x fliplr(x)], [y-CI fliplr(y+CI)], gray, 'FaceAlpha',0.2, 'EdgeColor','none')
clear p;

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;

a = gca;
if strcmp(wm, 'fa')
    a.YLabel.String = {'Fractional Anisotropy (FA)'};
elseif strcmp(wm, 'md')
    a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
end
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = {'Age'};
a.XLabel.FontSize = fontsize;
% title('Subregion: SUB')
title('')

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
pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['regression_' wm '_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;
 