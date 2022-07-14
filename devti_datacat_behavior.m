% Age - Microstructural Relationships Among HPC Subfields

clear all; close all; clc
format long g

% Set working directories.
indir = '/Volumes/Seagate/devti_devHPCsubfields/taskperf';
outdir = '/Volumes/Seagate/devti_devHPCsubfields/supportFiles';

%% Load behavioral data: task1.
temp = load(fullfile(indir, 'task1_behavior_n=78.mat'));

% Select data of interest and convert to table.
task1 = array2table(cat(2, temp.subnumbers', temp.age, temp.allPerfData, temp.allrtData), 'VariableNames', ...
    [{'subID'}, {'age'}, strcat('task1_acc_', temp.allPerfDataCols), strcat('task1_rt_', temp.allrtDataCols)]);
% Remove the extra group and age columns within the original task1 mat file
task1(:, [18 33 34]) = [];
% Move and rename the group column.
task1.Properties.VariableNames{17} = 'group';
task1 = [task1(:, 1:2) task1(:, 17) task1(:, 3:16) task1(:, 18:end)];

% Remove data identified as outliers by Schlichting et al., 2017. 
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
remove1 = unique([removetemp1 removetemp2]);

clear temp removetemp1 removetemp2;

% %% Load behavioral data: task3.
% temp = load(fullfile(indir, 'task3_behavior_n=78.mat'));
% 
% % Select data of interest and convert to table.
% task3 = array2table(cat(2, temp.includeSubsPerf', temp.meanacc), 'VariableNames', ...
%     [{'subID'}, {'task3_meanacc'}]);
% 
% % Remove data identified as outliers. 
% remove3 = temp.subnumbers(find(~ismember(temp.subnumbers, temp.includeSubsPerf)));
% 
% clear temp;
% 
%% Load behavioral data: task4.
temp = load(fullfile(indir, 'task4_behavior_n=78.mat'));

% Select data of interest and convert to table.
task4 = array2table(cat(2, temp.subnumbers', temp.rel_0acc, temp.rel_1acc, temp.rel_2acc, ...
    temp.rel_0rt, temp.rel_1rt, temp.rel_2rt), 'VariableNames', ...
    [{'subID'}, {'task4_rel0acc'}, {'task4_rel1acc'}, {'task4_rel2acc'}, {'task4_rel0rt'}, {'task4_rel1rt'}, {'task4_rel2rt'}]);

% Remove data identified as outliers. 
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
remove4 = unique([removetemp1 removetemp2]);

clear temp removetemp1 removetemp2;

%% Load behavioral data: task5.
temp = load(fullfile(indir, 'task5_behavior_n=78.mat'));

% Select data of interest and convert to table.
task5 = array2table(cat(2, temp.subnumbers', temp.GenDIRperf, temp.GenINDperf, temp.SpDIRperf, temp.SpINDperf, temp.ColDIRperf, temp.ColINDperf, ...
temp.GenDIRrt, temp.GenINDrt, temp.SpDIRrt, temp.SpINDrt, temp.ColDIRrt, temp.ColINDrt), 'VariableNames', ...
[{'subID'}, {'task5_GenDIRacc'}, {'task5_GenINDacc'}, {'task5_SpDIRacc'}, {'task5_SpINDacc'}, {'task5_ColDIRacc'}, {'task5_ColINDacc'}, ...
{'task5_GenDIRrt'}, {'task5_GenINDrt'}, {'task5_SpDIRrt'}, {'task5_SpINDrt'}, {'task5_ColDIRrt'}, {'task5_ColINDrt'}]);

% Remove data identified as outliers. 
removetemp1 = temp.subnumbers(find(temp.includeBoolPerf == 0));
removetemp2 = temp.subnumbers(find(temp.includeBoolRT == 0));
remove5 = unique([removetemp1 removetemp2]);

clear temp removetemp1 removetemp2;

% % Select subIDs to remove based on behavior. 
% remove = unique([remove1 remove5]);
% 
% keep = task1.subID(find(~ismember(task1.subID, remove)));
% task1 = task1(find(ismember(task1.subID, keep)), :); clear keep;
% 
% % keep = task3.subID(find(~ismember(task3.subID, remove)));
% % task3 = task3(find(ismember(task3.subID, keep)), :); clear keep;
% % 
% % keep = task4.subID(find(~ismember(task4.subID, remove)));
% % task4 = task4(find(ismember(task4.subID, keep)), :); clear keep;
% 
% keep = task5.subID(find(~ismember(task5.subID, remove)));
% task5 = task5(find(ismember(task5.subID, keep)), :); clear keep;

beh = [task1 task4(:, 2:end) task5(:, 2:end)];

% Save and export data.
filename = sprintf('devti_data_beh_%s', datestr(now,'yyyymmdd'));

% save it as a matlab table
save(fullfile(outdir, filename), 'beh')

% Save as a CSV files.
writetable(beh, fullfile(outdir, [filename '.csv']))
