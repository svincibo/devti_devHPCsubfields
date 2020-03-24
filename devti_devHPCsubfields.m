% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

blprojectid = '5e5672430f7fa65e1d3c9621';

% Select WM measure.
wm = {'fa', 'ad', 'rd', 'md'};
subregion = {'b_ca1', 'b_ca23', 'b_sub'};

% Load removals: statistical outliers.
load('devti_remove_statoutliers.mat')

% % Load removals: motion.
% load('devti_remove_motion.mat')
% 
% % Load removals: snr.
% load('devti_remove_snr.mat')

degp = 1;

% Make up the sex covariate for now.
% sex = Shuffle(cat(1, ones(29, 1), 2*ones(29, 1)));

for w = 1%:length(wm)
    
    load(['devti_data_' wm{w} '.mat'])
    
    for r = 1%:length(subregion)
        
        % Select outliers to remove.
        if strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_ca1')
            remove = outliers.fa_b_ca1;
        elseif strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_ca23')
            remove = outliers.fa_b_ca23;
        elseif strcmp(wm{w}, 'fa') && strcmp(subregion{r}, 'b_sub')
            remove = outliers.fa_b_sub;
        elseif strcmp(wm{w}, 'ad') && strcmp(subregion{r}, 'b_ca1')
            remove = outliers.ad_b_ca1;
        elseif strcmp(wm{w}, 'ad') && strcmp(subregion{r}, 'b_ca23')
            remove = outliers.ad_b_ca23;
        elseif strcmp(wm{w}, 'ad') && strcmp(subregion{r}, 'b_sub')
            remove = outliers.ad_b_sub;
        elseif strcmp(wm{w}, 'rd') && strcmp(subregion{r}, 'b_ca1')
            remove = outliers.rd_b_ca1;
        elseif strcmp(wm{w}, 'rd') && strcmp(subregion{r}, 'b_ca23')
            remove = outliers.rd_b_ca23;
        elseif strcmp(wm{w}, 'rd') && strcmp(subregion{r}, 'b_sub')
            remove = outliers.ad_b_sub;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_ca1')
            remove = outliers.md_b_ca1;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_ca23')
            remove = outliers.md_b_ca23;
        elseif strcmp(wm{w}, 'md') && strcmp(subregion{r}, 'b_sub')
            remove = outliers.md_b_sub;
        end
        
        % Mean center continuous variables.
        m = double(m - nanmean(m, 1));
        
        % Convert data to table for easier model specification.
        data = array2table(cat(2, transpose(sub), transpose(age), transpose(sex), transpose(iq), m), 'VariableNames', {'subID', 'age', 'sex',  'iq', roi{1, :}});
         
        % 1. Linear main effects model.
        modelspec = 'b_ca1 ~ sex + age + (1|subID)';
        if sum(remove) == 0          

            % Fit robust regression model.
            mdlr = fitlme(data, modelspec);
            
        else
            
            % Fit robust regression model, excluding outliers.
            mdlr = fitlme(data, modelspec, 'Exclude', remove);
            
        end
        
        % Correct AIC for sample size and predictor number: AICc.       
        aicc = mdlr.ModelCriterion.AIC + 2*size(mdlr.PredictorNames, 1)*((size(mdlr.PredictorNames, 1) + 1)/(size(mdlr.ObservationInfo, 1) - size(mdlr.PredictorNames, 1) - 1));
        
        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including linear main effects is ' num2str(aicc) '.']);
        
        % 2. Nonlinear main effects model.
        modelspec = 'b_ca1 ~ sex + age^2 + (1|subID)';
        if sum(remove) == 0          

            % Fit robust regression model.
            mdlr = fitlme(data, modelspec);
            
        else
            
            % Fit robust regression model, excluding outliers.
            mdlr = fitlme(data, modelspec, 'Exclude', remove);
            
        end
        
        % Correct AIC for sample size and predictor number: AICc.
        aicc = mdlr.ModelCriterion.AIC + 2*size(mdlr.PredictorNames, 1)*((size(mdlr.PredictorNames, 1) + 1)/(size(mdlr.ObservationInfo, 1) - size(mdlr.PredictorNames, 1) - 1));
        
        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including nonlinear main effects is ' num2str(aicc) '.']);
        
                % 3. Linear interaction model.
        modelspec = 'b_ca1 ~ sex*age + (1|subID)';
        if sum(remove) == 0          

            % Fit robust regression model.
            mdlr = fitlme(data, modelspec);
            
        else
            
            % Fit robust regression model, excluding outliers.
            mdlr = fitlme(data, modelspec, 'Exclude', remove);
            
        end
        
        % Correct AIC for sample size and predictor number: AICc.
        aicc = mdlr.ModelCriterion.AIC + 2*size(mdlr.PredictorNames, 1)*((size(mdlr.PredictorNames, 1) + 1)/(size(mdlr.ObservationInfo, 1) - size(mdlr.PredictorNames, 1) - 1));
        
        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including linear interactions is ' num2str(aicc) '.']);
        
        
        % 4. Nonlinear interactions model.
        modelspec = 'b_ca1 ~ sex*age + sex*(age^2) + (1|subID)';
        if sum(remove) == 0          

            % Fit robust regression model.
            mdlr = fitlme(data, modelspec);
            
        else
            
            % Fit robust regression model, excluding outliers.
            mdlr = fitlme(data, modelspec, 'Exclude', remove);
            
        end
        
        % Correct AIC for sample size and predictor number: AICc.
        aicc = mdlr.ModelCriterion.AIC + 2*size(mdlr.PredictorNames, 1)*((size(mdlr.PredictorNames, 1) + 1)/(size(mdlr.ObservationInfo, 1) - size(mdlr.PredictorNames, 1) - 1));
        
        % Get AIC value corrected for sample size: AICc.
        disp(['AICc for ' wm{w} ' in ' subregion{r} ' using a model including nonlinear interactions is ' num2str(aicc) '.']);
            
    end
    
    
    
    
    
    
    
    
    
    if strcmp(wm{w}, 'fa')
        
        ylab = 'Fractional Anisotropy';
        ymin = 0; ymax = 0.4;
        
    elseif strcmp(wm{w}, 'ad')
        
        ylab = {'Axial Diffusivity'; '(AD x 1e-3)'};
        ymin = 0.5; ymax = 1.5;
        
    elseif strcmp(wm{w}, 'rd')
        
        ylab = {'Radial Diffusivity'; '(RD x 1e-3)'};
        ymin = 0.5; ymax = 1;
        
    elseif strcmp(wm{w}, 'md')
        
        ylab = {'Mean Diffusivity'; '(MD x 1e-3)'};
        ymin = 0.5; ymax = 1;
        
    end
    
    %% Visualize.
    
    figure(w)
    
    % Correlations between age and average WM measure in ROI.
    % CA1
    x = transpose(age);
    y = m(strcmp(roi, 'b_ca1'));
    c1_color = [0 0.4470 0.7410];
    scatter(x, y, 'filled', 'MarkerEdgeColor', c1_color, 'MarkerFaceColor', c1_color)
    hold on;
    
    % c1 = corr(x(~isnan(y)), y(~isnan(y)));
    c1 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
    x1 = linspace(0,30);
    y1 = 1./(1+x1);
    
    f1 = polyval(c1,x1);
    % plot(x1,y1)
    plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c1_color)
    clear x
    hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
    hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c1_color(1:3));
    set(hp3, 'facecolor', c1_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
    
    disp(['CA1 = ' num2str(c1)])
    clear x y f1 hp3
    
    % CA23dg
    x = transpose(age);
    y = m(strcmp(roi, 'b_ca23dg'));
    c2_color = [0.4940 0.1840 0.5560];
    scatter(x(~isnan(y)), y(~isnan(y)), 'filled', 'MarkerEdgeColor', c2_color, 'MarkerFaceColor', c2_color)
    
    % c2 = corr(x(~isnan(y)), y(~isnan(y)));
    c2 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
    
    f1 = polyval(c2,x1);
    % plot(x1,y1)
    plot(x1,f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c2_color)
    clear x
    hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
    hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c2_color(1:3));
    set(hp3, 'facecolor', c2_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
    
    disp(['CA23dg = ' num2str(c2)])
    clear x y f1 hp3
    
    % SUB
    x = transpose(age);
    y = m(strcmp(roi, 'b_sub'));
    c3_color = [0.8500 0.3250 0.0980];
    scatter(x(~isnan(y)), y(~isnan(y)), 'filled', 'MarkerEdgeColor', c3_color, 'MarkerFaceColor', c3_color)
    
    % c3 = corr(x(~isnan(y)), y(~isnan(y)));
    c3 = polyfit(x(~isnan(y)), y(~isnan(y)), degp);
    
    f1 = polyval(c3,x1);
    % plot(x1,y1)
    plot(x1,f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', c3_color)
    clear x
    hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); x = (1:size(f1, 2))'.*.30;
    hp3 = patch([x; x(end:-1:1); x(1)], [lo'; hi(end:-1:1)'; lo(1)], c3_color(1:3));
    set(hp3, 'facecolor', c3_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
    
    disp(['SUB = ' num2str(c3)])
    
    legend([{'CA1'}, {''}, {''}, {'CA23dg'}, {''}, {''}, {'SUB'}, {''}], 'Location', 'southwest')
    legend box off
    
    ylabel(ylab)
    xlabel('Age (years)')
    
    xlim([0 30])
    ylim([ymin ymax])
    
    hold off;
    
    clear m roi sub x y f1 hp3
    
end

