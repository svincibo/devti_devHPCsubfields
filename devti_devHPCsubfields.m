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
subregion = {'b_ca1', 'b_ca23', 'b_dg', 'b_sub'};

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
       
    % Tell Matlab that sex and age group are categorical variables.
    data.sex = categorical(data.sex);
    data.group = categorical(data.group);
    
    % One subfield at a time.
    for r = 1:length(subregion)
        
        % Select outliers to remove.
        if strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_ca1')
            remove = outliers.fa_b_ca1;
        elseif strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_ca23')
            remove = outliers.fa_b_ca23;
        elseif strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_sub')
            remove = outliers.fa_b_sub;
        elseif strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_dg')
            remove = outliers.fa_b_dg;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_ca1')
            remove = outliers.md_b_ca1;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_ca23')
            remove = outliers.md_b_ca23;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_sub')
            remove = outliers.md_b_sub;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_dg')
            remove = outliers.md_b_dg;
        end
        
        % 1. Linear main effects model.
        modelspec = [subregion{r} ' ~ sex + age'];
        if sum(remove) == 0
            
            % Fit regression model.
            mdlr1 = fitlm(data, modelspec);
            
        else
            
            % Fit regression model, excluding outliers.
            mdlr1 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
            
        end

        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including linear main effects is ' num2str(mdlr1.ModelCriterion.AICc) '.']);
        disp(['Adjusted R-squared for ' wm{w} ' in ' subregion{r} ' using a model including linear main effects is ' num2str(mdlr1.Rsquared.Adjusted) '.']);
        
        mdlr1.anova
        clear mdlr
        
        % 2. Nonlinear main effects model.
        modelspec = [subregion{r} ' ~ sex + age^2'];
        if sum(remove) == 0
            
            % Fit regression model.
            mdlr2 = fitlm(data, modelspec);
            
        else
            
            % Fit regression model, excluding outliers.
            mdlr2 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
            
        end
 
        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including nonlinear main effects is ' num2str(mdlr2.ModelCriterion.AICc) '.']);
        disp(['Adjusted R-squared for ' wm{w} ' in ' subregion{r} ' using a model including nonlinear main effects is ' num2str(mdlr2.Rsquared.Adjusted) '.']);
        
        mdlr2.anova
        clear mdlr
        
        % 3. Linear interaction model.
        modelspec = [subregion{r} ' ~ sex*age'];
        if sum(remove) == 0
            
            % Fit regression model.
            mdlr3 = fitlm(data, modelspec);
            
        else
            
            % Fit regression model, excluding outliers.
            mdlr3 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
            
        end

        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including linear interactions is ' num2str(mdlr3.ModelCriterion.AICc) '.']);
        disp(['Adjusted R-squared for ' wm{w} ' in ' subregion{r} ' using a model including linear interactions is ' num2str(mdlr3.Rsquared.Adjusted) '.']);
        
        mdlr3.anova
        
        % 4. Nonlinear interactions model.
        modelspec = [subregion{r} ' ~ sex*(age^2)'];
        if sum(remove) == 0
            
            % Fit regression model.
            mdlr4 = fitlm(data, modelspec);
            
        else
            
            % Fit regression model, excluding outliers.
            mdlr4 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
            
        end

        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including nonlinear interactions is ' num2str(mdlr4.ModelCriterion.AICc) '.']);
        disp(['Adjusted R-squared for ' wm{w} ' in ' subregion{r} ' using a model including nonlinear interactions is ' num2str(mdlr4.Rsquared.Adjusted) '.']);
        
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
        
        % Known from prior testing: linear interactions model is best for FA
        % while non-linear interactions model is best for MD.
        if strcmp(wm{w}, 'fa')
            
            degp = 2;
            
            modelspec1 = [subregion{r} ' ~ sex'];
            modelspec2 = 'res ~ (age)^2';
            
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
            
            modelspec1 = [subregion{r} ' ~ sex'];
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
        
        % Set the color for this roi and add back in mean.
        if strcmp(subregion{r}, 'b_ca1')
            clr = [0 0.4470 0.7410]; % blue
%             y = y + nanmean(data_raw.b_ca1);
        elseif strcmp(subregion{r}, 'b_ca23')
            clr = [0.4940 0.1840 0.5560]; % purple
%             y = y + nanmean(data_raw.b_ca23);
        elseif strcmp(subregion{r}, 'b_dg')
            clr = [0.5019 0.5019 0.5019]; % gray
%             y = y + nanmean(data_raw.b_dg);
        elseif strcmp(subregion{r}, 'b_sub')
            clr = [0.8500 0.3250 0.0980]; % orange
%             y = y + nanmean(data_raw.b_sub);
        end
        
        if r == 1
            figure
            hold on;
        end
        
        % Correlations between age and average WM measure in ROI.
        scatter(x(~isnan(y) & sex == 1), y(~isnan(y) & sex == 1), 'Marker', 'x', 'MarkerEdgeColor', clr, 'LineWidth', linewidth)
        hold on;
        scatter(x(~isnan(y) & sex == 2), y(~isnan(y) & sex == 2), 'Marker', 'o', 'MarkerEdgeColor', clr, 'LineWidth', linewidth)
        
        c1 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
        x1 = linspace(0,30);
        
        f1 = polyval(c1,x1);
        plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', clr)
        hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x2 = (1:size(f1, 2))'.*.30;
        hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], clr);
        set(hp3, 'facecolor', clr, 'edgecolor', 'none', 'facealpha', .2);
        
        disp([subregion{r} ' = ' num2str(c1)])
        
    end
    
    legend([{'CA1, F'}, {'CA1, M'}, {''}, {''}, {'CA23, F'}, {'CA23, M'}, {''}, {''}, ...
        {'dg, M'}, {'dg, F'}, {''}, {''}, {'SUB, M'}, {'SUB, F'}, {''}, {''}], 'Location', 'bestoutside')
    legend box off
    
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
    
    print(fullfile(rootDir, 'plots', ['plot_linearfit_subfield_bygroup_'  wm{w}]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_linearfit_subfield_bygroup_' wm{w}]), '-depsc')
    
    hold off;
    
    clear m roi sub x y f1 hp3
    
    %% Perform a Two-way ANOVA for subfield: DV is wm measurement, Factor 1 is subfield: b_ca1, b_ca23, b_dg, b_sub, Factor 2 is age group: child, adolescent, adult
    
    % Transform into long form for glm using demeaned measurements.
    data_prefix = [data.subID data.age double(data.group) double(data.sex) data.iq];
    
    b_ca1 = data.b_ca1;
    b_ca1_idx = repmat(1, size(b_ca1));
    
    b_ca23 = data.b_ca23;
    b_ca23_idx = repmat(2, size(b_ca23));
    
    b_dg = data.b_dg;
    b_dg_idx = repmat(3, size(b_dg));
    
    b_sub = data.b_sub;
    b_sub_idx = repmat(4, size(b_sub));
    
    wm_measure = cat(1, b_ca1, b_ca23, b_dg, b_sub);
    wm_measure_idx = cat(1, b_ca1_idx, b_ca23_idx, b_dg_idx, b_sub_idx);
    data_prefix_new = repmat(data_prefix, [4 1]);
    
    d_array = cat(2, data_prefix_new, wm_measure_idx, wm_measure);
    d = array2table(d_array, 'VariableNames', {'subID', 'age', 'group', 'sex',  'iq', 'subfield', 'measurement'});
    
    % Tell Matlab that sex and age group are categorical variables.
    d.sex = categorical(d.sex);
    d.group = categorical(d.group);
    
    modelspec = 'measurement ~ sex + subfield*group';
    
    % Get outliers -- NOTE: for now this is any subject that was an outlier on any measure in any subfield.
    remove = sum(d.subID == unique(struct2array(outliers)), 2) >= 1;
    keep = sum(d.subID ~= unique(struct2array(outliers)), 2) >= 1;
    
    % Fit regression model, excluding outliers.
    mdlr = fitlm(d, modelspec, 'Exclude', remove);
    
    disp(wm{w})
    
    % Check for significant predictors.
    mdlr.anova
    
    clear data data_raw m_demeaned m
    
end

