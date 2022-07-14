% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

% Note bl-repository.
blprojectid = '5e5672430f7fa65e1d3c9621';

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Select measure and roi.
wm = 'md'; %'ad', 'rd', 'md'

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
head = [204 204 0]/255; % light burnt yellow
body = [204 190 0]/255; % burnt yellow
tail =  [210 43 43]/255; % cadmium red

hip = [2 129 129]/255; % dark turquoise

if strcmp(wm, 'fa')
    ylim_lo = 0.11; ylim_hi = 0.34;
elseif strcmp(wm, 'md')
    ylim_lo = .7; ylim_hi = 1.4;  
end


%% Hip

% Linear main effects model.
modelspec = 'tail ~ sex*age_cov^2';    

% Remove outliers.
removeidx = findoutliers(m, modelspec);
m2  = m; m2(removeidx, :) = []; 

% Mean center continuous variables for modelling.
m2(:, 5:end) = array2table(double(table2array(m2(:, 5:end)) - nanmean(table2array(m2(:, 5:end)), 1)));
mdl = fitlm(m2, modelspec) %, 'Exclude', find(sum(m.subID == remove, 2)));
n = mdl.NumObservations;
display(['N = ' num2str(n)])
display(['AICc = ' num2str(abs(mdl.ModelCriterion.AICc))])
clear m2;
 