% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

blprojectid = '5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields/';
yticklength = 0;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0.05;
linewidth = 2;
%%%%%%%%%%%%%%% TESTING AREA %%%%%%%%%%%%%%%%

% Use only the subjects that Meg used.
% data_meg = readtable(fullfile(rootDir, 'supportFiles/dti_subfields_long.csv'));
% sub_include = unique(data_meg.subnr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select WM measure.
wm = {'fa', 'md'}; %'ad', 'rd',
subregion = {'b_head'};

% Load removals: statistical outliers. (note motion outlier, sub90, is
% included in statoutliers)
load(fullfile(rootDir, 'supportFiles/devti_remove_statoutliers.mat'))

%
% % Load removals: motion.
% load('devti_remove_motionoutliers.mat')
%
% % Load removals: snr.
% load('devti_remove_snroutliers.mat')
% % Note: no datasets were removed due to unusually low snr

% One microstructural measurement at a time.
for w = 1:length(wm)
    
    load(fullfile(rootDir, ['supportFiles/devti_data_' wm{w} '_mrtrix3act.mat']))
    
    % Scale md, rd, and ad values for analysis and visualization.
    if ~strcmp(wm{w}, 'fa')
        m = m*1000;
    end
    
    %     % log transform
    %     m = log10(m);
    
    % Mean center continuous variables for modelling.
    m_demeaned = double(m - nanmean(m, 1));
    
    % Convert data to table for easier model specification.
    data = array2table(cat(2, transpose(sub), transpose(age), transpose(group), transpose(sex), transpose(iq), m_demeaned), 'VariableNames', {'subID', 'age', 'group', 'sex',  'iq', roi{1, :}});
    data_raw = array2table(cat(2, transpose(sub), transpose(age), transpose(group), transpose(sex), transpose(iq), m), 'VariableNames', {'subID', 'age', 'group', 'sex',  'iq', roi{1, :}});
    
    % Create columns for hippocampal head by averaging across all hippocampal subfields within the head.
    data.b_all = nanmean([data.b_ca1 data.b_ca23 data.b_dg data.b_sub], 2);
    data.l_all = nanmean([data.l_ca1 data.l_ca23 data.l_dg data.l_sub], 2);
    data.r_all = nanmean([data.r_ca1 data.r_ca23 data.r_dg data.r_sub], 2);
    
    % Tell Matlab that sex and age group are categorical variables.
    data.sex = categorical(data.sex);
    data.group = categorical(data.group);

    % Select outliers to remove.
    if strcmp(wm{w}, 'fa')
        remove = unique(cat(2, outliers.fa_b_ca1, outliers.fa_b_ca23, outliers.fa_b_sub, outliers.fa_b_dg));
    elseif strcmp(wm{w}, 'md')
        remove = unique(cat(2, outliers.md_b_ca1, outliers.md_b_ca23, outliers.md_b_sub, outliers.md_b_dg));
    end
    
    % 1. Linear main effects model.
    modelspec = 'b_all ~ sex + age';
    if sum(remove) == 0
        
        % Fit regression model.
        mdlr1 = fitlm(data, modelspec);
        
    else
        
        % Fit regression model, excluding outliers.
        mdlr1 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
        
    end

    % Get AIC value corrected for sample size: AICc.
    disp(['AICc for ' wm{w} ' in hippocampal head using a model including linear main effects is ' num2str(mdlr1.ModelCriterion.AICc) '.']);
    disp(['Adjusted R-squared for ' wm{w} ' in hippocampal head using a model including linear main effects is ' num2str(mdlr1.Rsquared.Adjusted) '.']);
    
    mdlr1.anova
    clear mdlr
    
    % 2. Nonlinear main effects model.
    modelspec = 'b_all ~ sex + age^2';
    if sum(remove) == 0
        
        % Fit regression model.
        mdlr2 = fitlm(data, modelspec);
        
    else
        
        % Fit regression model, excluding outliers.
        mdlr2 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
        
    end

    % Get AIC value corrected for sample size: AICc.
    disp(['AICc for ' wm{w} ' in hippocampal head using a model including nonlinear main effects is ' num2str(mdlr2.ModelCriterion.AICc) '.']);
    disp(['Adjusted R-squared for ' wm{w} ' in hippocampal head using a model including nonlinear main effects is ' num2str(mdlr2.Rsquared.Adjusted) '.']);
    
    mdlr2.anova
    clear mdlr
    
    % 3. Linear interaction model.
    modelspec = 'b_all ~ sex*age';
    if sum(remove) == 0
        
        % Fit regression model.
        mdlr3 = fitlm(data, modelspec);
        
    else
        
        % Fit regression model, excluding outliers.
        mdlr3 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
        
    end

    % Get AIC value corrected for sample size: AICc.
    disp(['AICc for ' wm{w} ' in hippocampal head using a model including linear interactions is ' num2str(mdlr3.ModelCriterion.AICc) '.']);
    disp(['Adjusted R-squared for ' wm{w} ' in hippocampal head using a model including linear interactions is ' num2str(mdlr3.Rsquared.Adjusted) '.']);
    
    mdlr3.anova
    
    % 4. Nonlinear interactions model.
    modelspec = 'b_all ~ sex*(age^2)';
    if sum(remove) == 0
        
        % Fit regression model.
        mdlr4 = fitlm(data, modelspec);
        
    else
        
        % Fit regression model, excluding outliers.
        mdlr4 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
        
    end
  
    % Get AIC value corrected for sample size: AICc.
    disp(['AICc for ' wm{w} ' in hippocampal head using a model including nonlinear interactions is ' num2str(mdlr4.ModelCriterion.AICc) '.']);
    disp(['Adjusted R-squared for ' wm{w} ' in hippocampal head using a model including nonlinear interactions is ' num2str(mdlr4.Rsquared.Adjusted) '.']);
    
    mdlr4.anova
    
    %% Visualize.
    
    %     % Subselect for testing: uncomment this when using only Meg's data points.
    %     if exist('sub_include')
    %         keep_idx = ismember(sub', sub_include);
    %     else
    %         keep_idx = 1:length(sub);
    %     end
    
    % Remove outliers so that they are not in the plot.
    if sum(remove) == 0
        keep_idx = 1:length(data.subID);
    else
        keep_idx = ~ismember(data.subID, remove);
    end
    
    % Select x data.
    x = data.age(keep_idx);
    
    % Select covariates.
    sex = double(data.sex(keep_idx));
    
    % Known from prior testing: linear interactions model (mdlr3) is best for FA
    % while non-linear interactions model is best for MD.
    if strcmp(wm{w}, 'fa')
        
        degp = 2;
        
        modelspec1 = 'b_all ~ sex + age';
        modelspec2 = 'res ~ age^2';
        
        if sum(remove) == 0
            
            % Get and remove residuals.
            mdlr = fitlm(data, modelspec1);
            data.res = table2array(mdlr.Residuals(:, 1));
            
            % Fit regression model.
            mdlr3 = fitlm(data, modelspec2);
            
        else
            
            data = data(keep_idx, :);
            
            % Get and remove residuals.
            mdlr = fitlm(data, modelspec1);
            data.res = table2array(mdlr.Residuals(:, 1));
            
            % Fit regression model, excluding outliers.
            mdlr3 = fitlm(data, modelspec2);
            
        end
        
        % Select y-data: ".. adjusted y-data were calculated by adding the residual
        % to the fitted value for each point." Figure 2 in Schlichting et al., JoCN, 2018.
        % I used the "raw residuals", i.e., column 1.
        y = mdlr3.Fitted + table2array(mdlr3.Residuals(:, 1));
        
    elseif strcmp(wm{w}, 'md')
        
        degp = 1;
        
        modelspec1 = 'b_all ~ sex';
        modelspec2 = 'res ~ age';
        
        if sum(remove) == 0
            
            % Get and remove residuals.
            mdlr = fitlm(data, modelspec1);
            data.res = table2array(mdlr.Residuals(:, 1));
            
            % Fit regression model.
            mdlr3 = fitlm(data, modelspec2);
            
        else
            
            data = data(keep_idx, :);
            
            % Get and remove residuals.
            mdlr = fitlm(data, modelspec1);
            data.res = table2array(mdlr.Residuals(:, 1));
            
            % Fit regression model, excluding outliers.
            mdlr3 = fitlm(data, modelspec2);
            
        end
        
        % Select y-data: ".. adjusted y-data were calculated by adding the residual
        % to the fitted value for each point." Figure 2 in Schlichting et al., JoCN, 2018.
        % I used the "raw residuals", i.e., column 1.
        y = mdlr3.Fitted + table2array(mdlr3.Residuals(:, 1));
        
    else
        
        degp = 1;
        
    end
    
    %         % Un-demean for visualization.
    %         y = y + double(m - nanmean(m, 1));
    
    % Set the color.
    clr = [0 0 0]; % black
    
    figure
    hold on;
    
    plotAdjustedResponse(mdlr3,'age')
% plotAdded(mdlr3, 3)
    
%     % Correlations between age and average WM measure in ROI.
%     scatter(x(~isnan(y) & sex == 1), y(~isnan(y) & sex == 1), 'Marker', 'x', 'MarkerEdgeColor', clr, 'LineWidth', linewidth)
%     hold on;
%     scatter(x(~isnan(y) & sex == 2), y(~isnan(y) & sex == 2), 'Marker', 'o', 'MarkerEdgeColor', clr, 'LineWidth', linewidth)
%     
%     % Plot fit for all age groups.
%     c1 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
%     x1 = linspace(0,30);
%     
%     f1 = polyval(c1,x1);
%     plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', clr)
%     hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x2 = (1:size(f1, 2))'.*.30;
%     hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], clr);
%     set(hp3, 'facecolor', clr, 'edgecolor', 'none', 'facealpha', .2);
%     
%     disp(['all subfields = ' num2str(c1)])
%     
%     legend([{'all subfields, M'}, {'all subfields, F'}, {'fit'}, {'sd'}], 'Location', 'southeast')
%     legend box off
    
    if strcmp(wm{w}, 'fa')
        
        ylab = 'Fractional Anisotropy, demeaned, adjusted';
        ylim_lo = -0.05; ylim_hi = 0.05;
        %         ylim_lo = -01; ylim_hi = -0.4; % for log transform
        
    elseif strcmp(wm{w}, 'md')
        
        ylab = {'Mean Diffusivity, (demeaned, adjusted, MD x 1e-3)'};
        ylim_lo = -0.1; ylim_hi = 0.1;
        %         ylim_lo = -0.2; ylim_hi = 0.2; % for log transform
        
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
    xax.FontAngle = fontangle;
    
    % yaxis
    yax = get(gca,'yaxis');
    yax.Limits = [ylim_lo ylim_hi];
    yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
    yax.TickDirection = 'out';
    yax.TickLength = [xticklength xticklength];
    yax.TickLabels = {num2str(ylim_lo, '%1.2f'), num2str((ylim_lo+ylim_hi)/2, '%1.0f'), num2str(ylim_hi, '%1.2f')};
    yax.FontName = fontname;
    yax.FontSize = fontsize;
    
    % general
    a = gca;
    %     set(gca, 'YScale', 'log')
    
    %     a.TitleFontWeight = 'normal';
    box off
    
    a.YLabel.String = ylab;
    a.YLabel.FontSize = fontsize;
    a.XLabel.String = 'Age (years)';
    
    pbaspect([1 1 1])
    
    print(fullfile(rootDir, 'plots', ['plot_linearfit_allsubfieldscombined_'  wm{w}]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_linearfit_allsubfieldscombined_' wm{w}]), '-depsc')
    
    hold off;
    
    
    %% Perform a One-way ANOVA for subfield: DV is wm measurement in head. Factor is age group: child, adolescent, adult
       
    modelspec = 'b_all ~ sex*group';
    
    % Get outliers -- NOTE: for now this is any subject that was an outlier on any measure in any subfield.
    remove = sum(data.subID == unique(struct2array(outliers)), 2) >= 1;
    keep = sum(data.subID ~= unique(struct2array(outliers)), 2) >= 1;
    
    % Fit regression model, excluding outliers.
    mdlr = fitlm(data, modelspec, 'Exclude', remove);
    
    disp(wm{w})
    
    % Check for significant predictors.
    mdlr.anova
    
    clear m roi sub x y f1 hp3
    
end