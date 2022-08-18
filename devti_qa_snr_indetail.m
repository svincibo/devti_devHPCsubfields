% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

remove_outliers = 'yes';
include = 'all'; % all, childrenonly
outlier = [11 90]; 

% Read in behavioral data.
behdata = readtable(fullfile(rootdir, 'supportFiles', 'devti_data_beh_forSPSS_20220705.csv'), 'TreatAsEmpty', {'.', 'na'});
data = table2array(behdata); clear behdata;

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behavioral/demographic data.
    if ~isempty(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        if ~isempty(sub_contents_snr)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get standard deviation of SNR in X, Y, and Z directions.
            sd(sub_count) = std([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));

            % Get ylabel.
            lab{sub_count} = grp_contents(s).name(5:7);

            % Get training group.
            group(sub_count) = data(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))), 3);

            % Get age.
            age(sub_count) = data(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))), 2);

            clear data_snr_temp get_temp

        end % if ~isempty

    end % end if exist
    
end % end s

% Concatenate and sort according to group.
toplot = cat(2, subID', group', m', sd', b0', age');
[toplot idxsort] = sortrows(toplot, [2 3 1], 'ascend');
lab = lab(:, idxsort);

% cheap trick
clear subID group m sd b0 age
subID = toplot(:, 1)';
group = toplot(:, 2)';
m = toplot(:, 3)';
sd = toplot(:, 4)';
b0 = toplot(:, 5)';

% SNR plot
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

yc_color  = [50 180 100]/255; 
oc_color = [50 100 180]/255; 
a_color = [100 50 180]/255;

% yc_color = [0.6350 0.0780 0.1840]; %red
% oc_color = [0 0.4470 0.7410]; %blue
% a_color = [0.41176 0.41176 0.41176]; %gray

figure(1); hold on;

gscatter(m, 1:length(m), group, [yc_color; oc_color; a_color], '.', 20)
gscatter(b0, 1:length(b0), group, [yc_color; oc_color; a_color], 'x', 8)
% plot([5 5], [0 length(subID)+0.5], ':k')

for p = 1:length(m)
    
    if group(p) == 1
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', yc_color)
        
    elseif group(p) == 2
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', oc_color)
        
    elseif group(p) == 3
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', a_color)
              
    end
    
end
gscatter(m, 1:length(m), group, [yc_color; oc_color; a_color], '.', 20)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [min([b0 m + sd])-5 max([b0 m + sd])+5];
xax.TickValues = floor(min([b0 m + sd])):10:ceil(max([b0 m + sd]));
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(subID)+0.5];
yax.TickValues = 1:length(subID);
yax.TickDirection = 'out';
yax.TickLabels = lab;
yax.FontName = fontname;
yax.FontSize = 6;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adolescents', 'Adults'}, 'Location', 'northeast');
legend box off

a.XLabel.String = 'SNR';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])

a.YLabel.String = 'Subject ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';

print(fullfile(rootdir, 'plots', 'plot_snr'), '-dpng')

hold off;

%% plot correlation with age

toplot = array2table(toplot);
toplot.Properties.VariableNames = {'subID', 'age_group', 'snr_dwi', 'snr_dwi_sd', 'snr_b0', 'age_cov'};

% Linear main effects model.
modelspec = 'snr_dwi ~ age_cov';  

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
a.YLabel.String = {'SNR, weighted volumes'};
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
ylim_lo = 0; ylim_hi = 45;
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [5 15 25 35 45];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
xlabels = {'5', '15', '25', '35', '45'};
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
text(6, 44, ['adjr2 = ' num2str(mdl.Rsquared.Ordinary)])
text(6, 42, ['rmse = ' num2str(mdl.RMSE)])
text(6, 40, ['beta = ' num2str(mdl.Coefficients.Estimate(2))])
text(6, 38, ['p = ' num2str(mdl.Coefficients.pValue(2))])
text(6, 36, ['n = ' num2str(n)])

print(fullfile(rootdir, 'plots', ['qa_snr_' modelspec '_n=' num2str(n)]), '-dpng')
print(fullfile(rootdir, 'plots', 'eps', ['qa_snr_' modelspec '_n=' num2str(n)]), '-depsc')

hold off;

% Save snr.
filename = sprintf('devti_qa_snr_%s', datestr(now,'yyyymmdd'));

% save it as a matlab table
save(fullfile(rootdir, 'supportFiles', filename), 'toplot')

% Save as a CSV files.
writetable(toplot, fullfile(rootdir, 'supportFiles', [filename '.csv']))