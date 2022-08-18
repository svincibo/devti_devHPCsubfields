% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

% Note bl-repository.
blprojectid = '5e5672430f7fa65e1d3c9621';

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Load snr data.
datestring = '20220810';
filename = sprintf('devti_qa_snr_%s', datestring);
load(fullfile(rootdir, 'supportFiles', filename))
snr = toplot; clear toplot;

% Load snr data.
datestring = '20220810';
filename = sprintf('devti_qa_fd_%s', datestring);
load(fullfile(rootdir, 'supportFiles', filename))
fd = toplot; clear toplot;

% Load the white matter data.
wm = 'fa'; %'ad', 'rd', 'md'
datestring = '20220705';
filename = sprintf('devti_data_%s_%s', wm, datestring);
load(fullfile(rootdir, 'supportFiles', filename))
mri = m; clear m;

% Rectify subIDs.
idx = find(ismember(mri.subID, fd.subID));
mri = mri(idx, :); clear idx;
idx = find(ismember(fd.subID, mri.subID));
fd = fd(idx, :); clear idx;
idx = find(ismember(snr.subID, mri.subID));
snr = snr(idx, :); clear idx;

% Concatenate beh and mri data into one table.
idx_fd = find(strcmp(fd.Properties.VariableNames, 'fd'));
idx_snr = find(strcmp(snr.Properties.VariableNames, 'snr_dwi'));
m = [mri fd(:, idx_fd) snr(:, idx_snr)];

% Scale md values for analysis and visualization, if this is md.
if strcmp(wm, 'md')
    scalefactor = 1000;
elseif strcmp(wm, 'fa')
    scalefactor = 1;
end
m(:, 7:22) = array2table(table2array(m(:, 7:22)).*scalefactor);

% % Sort according to age.
% age_idx = find(strcmp(m.Properties.VariableNames, 'age_cov') == 1);
% m = sortrows(m, age_idx);

% Select rois.
rois= m.Properties.VariableNames(7:22);

% Indicate qa measures to target.
qa= {'fd', 'snr'};

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

count = 0;
for q = 1:length(qa)

    for r = 1:length(rois)

        count = count + 1;

        % Find the column that contains the data for the current roi and extract the data.
        idx_roi = find(strcmp(m.Properties.VariableNames, rois{r}) == 1);
        roidata = table2array(m(:, idx_roi));

        if strcmp(qa{q}, 'fd')
            toplot = array2table([m.subID m.fd m.age_group roidata]);
        elseif strcmp(qa{q}, 'snr')
            toplot = array2table([m.subID m.snr_dwi m.age_group roidata]);
        end
        toplot.Properties.VariableNames = {'subID', 'measure', 'age_group', 'roidata'};

        figure(count); hold on;
        set(gcf,'Visible','on');
        set(0,'DefaultFigureVisible','on');

        % Linear main effects model.
        modelspec = 'roidata ~ measure';

        % % Remove outliers, use only *the data that were used*. Were removed before
        % % this step.
        % removeidx = find(sum(toplot.subID == outlier, 2));
        toplot2 = toplot; %toplot2(removeidx, :) = [];

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
        toplot2(:, 4) = array2table(double(table2array(toplot2(:, 4)) - nanmean(table2array(toplot2(:, 4)), 1)));
        mdl = fitlm(toplot2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
        n = mdl.NumObservations;
        display(['N = ' num2str(n)])
        display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
        clear toplot2;

        a = gca;
        if strcmp(wm, 'md')
            a.YLabel.String = {'Mean Diffusivity (MD)'; '(MD x 1e-3)'};
        elseif strcmp(wm, 'fa')
            a.YLabel.String = {'Fractional Anisotropy (FA)'};
        end
        a.YLabel.FontSize = fontsize;
        a.YLabel.FontAngle = fontangle;
        if strcmp(qa{q}, 'fd')
            a.XLabel.String = {'Framewise Displacement (FD)'};
        elseif strcmp(qa{q}, 'snr')
            a.XLabel.String = {'SNR, weighted volumes'};
        end
        a.XLabel.FontSize = fontsize;
        % title('Subregion: SUB')
        title(rois{r})

        % xaxis
        if strcmp(qa{q}, 'fd')
            xlim_lo = 0.40; xlim_hi = 1.00;
        elseif strcmp(qa{q}, 'snr')
            xlim_lo = 0; xlim_hi = 25;
        end
        xax = get(gca, 'xaxis');
        xax.Limits = [xlim_lo xlim_hi];
        xax.TickValues = [xlim_lo (xlim_lo+xlim_hi)/2 xlim_hi];
        xax.TickDirection = 'out';
        xax.TickLength = [yticklength yticklength];
        if strcmp(qa{q}, 'fd')
            xlabels = {num2str(xlim_lo, '%1.1f'), num2str((xlim_lo+xlim_hi)/2, '%1.1f'), num2str(xlim_hi, '%1.1f')};
        elseif strcmp(qa{q}, 'snr')
            xlabels = {num2str(xlim_lo, '%1.0f'), num2str((xlim_lo+xlim_hi)/2, '%1.0f'), num2str(xlim_hi, '%1.0f')};
        end
        xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
        xax.TickLabels = xlabels;
        xax.FontName = fontname;
        xax.FontSize = fontsize;

        % yaxis
        if strcmp(wm, 'md')
            ylim_lo = 0.6; ylim_hi = 1.3;
        elseif strcmp(wm, 'fa')
            ylim_lo = 0.11; ylim_hi = 0.34;
        end
        yax = get(gca,'yaxis');
        yax.Limits = [ylim_lo ylim_hi];
        yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
        yax.TickDirection = 'out';
        yax.TickLength = [xticklength xticklength];
        yax.TickLabels = {num2str(ylim_lo, '%1.1f'), num2str((ylim_lo+ylim_hi)/2, '%1.1f'), num2str(ylim_hi, '%1.1f')};
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
        if strcmp(wm, 'fa')
            if strcmp(qa{q}, 'snr')
                text(1, ylim_lo+.006, ['adjr2 = ' num2str(mdl.Rsquared.Ordinary)])
                text(1, ylim_lo+.016, ['rmse = ' num2str(mdl.RMSE)])
                text(1, ylim_lo+.026, ['beta = ' num2str(mdl.Coefficients.Estimate(2))])
                text(1, ylim_lo+.036, ['p = ' num2str(mdl.Coefficients.pValue(2))])
                text(1, ylim_lo+.046, ['n = ' num2str(n)])
            elseif strcmp(qa{q}, 'fd')
                text(0.42, ylim_hi-.01, ['adjr2 = ' num2str(mdl.Rsquared.Ordinary)])
                text(0.42, ylim_hi-.02, ['rmse = ' num2str(mdl.RMSE)])
                text(0.42, ylim_hi-.03, ['beta = ' num2str(mdl.Coefficients.Estimate(2))])
                text(0.42, ylim_hi-.04, ['p = ' num2str(mdl.Coefficients.pValue(2))])
                text(0.42, ylim_hi-.05, ['n = ' num2str(n)])
            end
        elseif strcmp(wm, 'md')
            if strcmp(qa{q}, 'snr')
                text(1, ylim_lo+.05, ['adjr2 = ' num2str(mdl.Rsquared.Ordinary)])
                text(1, ylim_lo+.07, ['rmse = ' num2str(mdl.RMSE)])
                text(1, ylim_lo+.09, ['beta = ' num2str(mdl.Coefficients.Estimate(2))])
                text(1, ylim_lo+.11, ['p = ' num2str(mdl.Coefficients.pValue(2))])
                text(1, ylim_lo+.13, ['n = ' num2str(n)])
            elseif strcmp(qa{q}, 'fd')
                text(0.42, ylim_hi-.05, ['adjr2 = ' num2str(mdl.Rsquared.Ordinary)])
                text(0.42, ylim_hi-.07, ['rmse = ' num2str(mdl.RMSE)])
                text(0.42, ylim_hi-.09, ['beta = ' num2str(mdl.Coefficients.Estimate(2))])
                text(0.42, ylim_hi-.11, ['p = ' num2str(mdl.Coefficients.pValue(2))])
                text(0.42, ylim_hi-.13, ['n = ' num2str(n)])
            end
        end

        print(fullfile(rootdir, 'plots', ['qa_' qa{q} '_' wm '_' rois{r} '_' modelspec '_n=' num2str(n)]), '-dpng')
        print(fullfile(rootdir, 'plots', 'eps', ['qa_' qa{q} '_' wm '_' rois{r} '_' modelspec '_n=' num2str(n)]), '-depsc')

        hold off;

    end % end roi

end % end qa
