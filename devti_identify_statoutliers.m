% "Remove multivariate outliers for each subregion by applying the box
% plot rule (Frigge, Hoaglin, & Iglewicz, 1989) to observations with
% unusually low weights computed using robust regression." (Schlichting
% et al., 2018)

clear all; close all; clc
format long g

load('/Volumes/240/devti_devHPCsubfields/supportFiles/devti_data_md_mrtrix3act.mat')

% Independent Variable
y = m(strcmp(roi, 'b_dg'));
        
% Dependent Variables
x = transpose(age);

% Fit robust regression model.
mdlr = fitlm(x, y, 'RobustOpts', 'on');
        
% % Examine model residuals: histogram of raw residuals.
% figure(1); plotResiduals(mdlr)

% Examine model residuals: boxplot of raw residuals.
figure(1); boxplot(mdlr.Residuals.Raw)
ylabel('Raw Residuals')

% Examine robust weights: boxplot of robust weights.
figure(2); boxplot(mdlr.Robust.Weights)
ylabel('Robust Beta-Weights')

% % Examine the model residuals: normal probability plot of raw residuals.
% figure(3); plotResiduals(mdlr, 'probability')
        
%         % Find index of outlier.
%         [~, outlier] = min(mdlr.Residuals.Raw); mdlr.Robust.Weights(outlier)
%         
%         % Check the median weight.
%         median(mdlr.Robust.Weights)
               
% % Manually record outliers. Include observations with outlier residuals
% and unusually low weights computed using robust regression. (0 indicates
% no outliers). 90 is included because it is an outlier in motion.
% Also: subs 6 and 90 had FA <=1.15 which is too low.
outliers.fa_b_ca1 = [6 90];
outliers.fa_b_ca23 = [6 16 51 90];
outliers.fa_b_sub = [6 13 90];
outliers.fa_b_dg = [6 70 74 90];

outliers.ad_b_ca1 = [6 90];
outliers.ad_b_ca23 = [6 90];
outliers.ad_b_sub = [6 90];
outliers.ad_b_dg = [6 90];

outliers.rd_b_ca1 = [6 90];
outliers.rd_b_ca23 = [6 90];
outliers.rd_b_sub = [6 90];
outliers.rd_b_dg = [6 90];

outliers.md_b_ca1 = [6 90];
outliers.md_b_ca23 = [6 90];
outliers.md_b_sub = [6 90];
outliers.md_b_dg = [6 90];

% save('/Volumes/240/devti_devHPCsubfields/supportFiles/devti_remove_statoutliers.mat', 'outliers')


    
