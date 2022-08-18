clear all; close all; clc
format shortG

yc_color  = [50 180 100]/255; 
oc_color = [50 100 180]/255; 
a_color = [100 50 180]/255;

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';
% addpath(genpath(fullfile(rootDir, 'proj-5e5672430f7fa65e1d3c9621')));

remove_outliers = 'yes';
include = 'all'; % all, childrenonly
outlier = [11 90]; 

%%%%%%%%%%%%%%% TESTING AREA %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in behavioral data.
beh_data_in_tbl = readtable([rootdir 'supportFiles/devti_data_beh_forSPSS_20220705.csv'], 'TreatAsEmpty', {'.', 'na'});

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
            % Remove the '.' and '..' files.
            sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
            
            % Get motion parameters for this subject.
            load(fullfile(sub_contents_motion.folder,sub_contents_motion.name));
            
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
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name(5:7);
            
            % Get training group.
            group(sub_count) = beh_data_in_tbl.group(s);  

            % Get age.
            age(sub_count) = beh_data_in_tbl.age(s);  
    
end % end s

meanmotion = mean(motion, 2);
sd = std(motion, 0, 2);

% Concatenate and sort according to group.
toplot = cat(2, subID', group', meanmotion, sd, age');
[toplot idxsort] = sortrows(toplot, [2 3 1], 'ascend');
lab = lab(:, idxsort);

% FD plot
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
xticklength = 0;
alphablend = .8;
yticklength = 0;
xticklength = 0.02;

figure(1); hold on; 
gscatter(toplot(:, 3), 1:length(toplot(:, 3)), toplot(:, 2), [yc_color; oc_color; a_color], '.', 30)

for p = 1:length(toplot(:, 3))

    if toplot(p, 2) == 1
        color = yc_color;
    elseif toplot(p, 2) == 2
        color = oc_color;
    elseif toplot(p, 2) == 3
        color = a_color;
    end
    plot([toplot(p, 3) - abs(toplot(p, 4)) toplot(p, 3) + abs(toplot(p, 4))], [p p], 'Color', color)

end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 2];
xax.TickValues = 0:0.5:2;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(lab)+0.5];
yax.TickValues = 1.5:2:ceil(length(lab));
yax.TickDirection = 'out';
yax.TickLabels = lab;
yax.FontName = fontname;
yax.FontSize = 8;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adolescents', 'Adults'}, 'Location', 'northeast');
legend box off

a.XLabel.String = 'Framewise Displacement (FD)';
a.XLabel.FontSize = fontsize;

a.YLabel.String = 'Subject ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';

pbaspect([1 1 1])

print(fullfile(rootdir, 'plots', 'plot_fd'), '-dpng')

hold off;

%% plot correlation with age

toplot = array2table(toplot);
toplot.Properties.VariableNames = {'subID', 'age_group', 'fd', 'fd_sd', 'age_cov'};

% Linear main effects model.
modelspec = 'fd ~ age_cov';  

% Remove outliers, need to show that there was no correlation between fd and age in *the data that were used*.
removeidx = find(sum(toplot.subID == outlier, 2));
toplot2  = toplot; toplot2(removeidx, :) = []; 
    
figure(2);
% Fit regression model.
mdl = fitlm(toplot2, modelspec); %, 'Exclude', find(sum(m.subID == remove, 2)));
p = plot(mdl);
hold on;

color = [128 128 128]/255;
p(1).Color = color;
p(1).Marker = 'o';
p(1).MarkerFaceColor = 'w';
p(1).MarkerEdgeColor = 'w';
p(1).MarkerSize = markersize;
p(2).Color = color;
p(2).LineWidth = linewidth;
p(3).LineStyle = 'none';
p(4).LineStyle = 'none';
% Plot data with group color coded.
s = gscatter(p(1).XData, p(1).YData, toplot2.age_group, [yc_color; oc_color; a_color]);
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
pci(1).Marker = 'none'; pci(2).Color = color;
pci(3).LineStyle = 'none'; pci(4).LineStyle = 'none'; 
clear p s pci;

% Mean center continuous variables for modelling.
toplot2(:, 5:end) = array2table(double(table2array(toplot2(:, 5:end)) - nanmean(table2array(toplot2(:, 5:end)), 1)));
mdl = fitlm(toplot2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear toplot2;

a = gca;
a.YLabel.String = {'Framewise Displacement (FD)'};
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
ylim_lo = 0; ylim_hi = 2;
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

% text(6, 1.8, ['aic/aicc = ' num2str(mdl.ModelCriterion.AIC) ', ' num2str(mdl.ModelCriterion.AICc)])
% text(6, 1.7, ['bic = ' num2str(mdl.ModelCriterion.BIC)])
text(6, 1.8, ['adjr2 = ' num2str(mdl.Rsquared.Ordinary)])
text(6, 1.7, ['rmse = ' num2str(mdl.RMSE)])
text(6, 1.6, ['beta = ' num2str(mdl.Coefficients.Estimate(2))])
text(6, 1.5, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(6, 1.4, ['n = ' num2str(n)])

print(fullfile(rootdir, 'plots', ['qa_fd_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['qa_fd_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

% Save fd.
filename = sprintf('devti_qa_fd_%s', datestr(now,'yyyymmdd'));

% save it as a matlab table
save(fullfile(rootdir, 'supportFiles', filename), 'toplot')

% Save as a CSV files.
writetable(toplot, fullfile(rootdir, 'supportFiles', [filename '.csv']))