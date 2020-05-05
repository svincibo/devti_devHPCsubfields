% This script reads in the tract profiles from Brad Caron's Tract Profile
% App and makes them into a group-level tract-profile.

clear all; close all; clc
format long g

blprojectid = 'proj-5e5672430f7fa65e1d3c9621/act';
% For fsldtifit measurements use: proj-5e5672430f7fa65e1d3c9621
% For mrtrix3 act measurements use: proj-5e5672430f7fa65e1d3c9621/act

% Select WM measure.
wm = {'fa', 'md'};

% Set working directories.
rootDir = '/Volumes/240/devti_devHPCsubfields';
%rng default % for reproducibility

% Get age.
load(fullfile(rootDir, 'supportFiles/data.mat'))
subs = array2table(data, 'VariableNames', {'subID', 'cov_age', 'iq', 'gp_age', 'a', 'b', 'c', 'd', 'e'});

% Add fix so to account for leading zeros in file names.
for k = 1:length(subs.subID)
    
    if subs.subID(k) < 10
        
        subID{k, 1} = strcat('00', num2str(subs.subID(k)));
        
    else
        
        subID{k, 1} = strcat('0', num2str(subs.subID(k)));
        
    end
    
end

%% ROI.

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Do for each wm measure.

for w = 1:size(wm, 2)
    
    % Load in each ROI's measures for this subject.
    sub_count = 0;
    for i = 1:size(grp_contents, 1)
        
        % Only read in data for subjects that we want to consider in this analysis.
        if sum(ismember(subID, grp_contents(i).name(end-2:end))) ~= 0
            
            disp(grp_contents(i).name(end-2:end))
            
            sub_count = sub_count +1;
            
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_rois = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*.nii.gz']));
            
            % Remove the '.' and '..' files.
            sub_contents_rois = sub_contents_rois(arrayfun(@(x) x.name(1), sub_contents_rois) ~= '.');
            
            for j = 1:size(sub_contents_rois)
                
                % Read in data for this subject and this tract.
                data_temp = niftiread(fullfile(sub_contents_rois(j).folder, sub_contents_rois(j).name));
                
                % Convert all zeros to NaN.
                data_temp(data_temp == 0) = NaN;
                
                % Average and standard deviation all non-NaN values.
                m(sub_count, j) = nanmean(nanmean(nanmean(data_temp, 1), 2), 3);
                %sd{sub_count, j} = nanstd(nanstd(nanstd(data_temp, 1), 2), 3);
                
                % Grab roi name.
                roi{sub_count, j} = sub_contents_rois(j).name(12:end-7);
                
                % Grab subID.
                sub(sub_count) = str2num(grp_contents(i).name(end-2:end));
                
                % Grab age.
                age(sub_count) = subs.cov_age(subs.subID == sub(sub_count));
                
                % Grab sex. 1 = F, 2 = M
                %sex(sub_count) = subs.cov_sex(subs.subID == sub(sub_count));
                sex(sub_count) = randi([1 2]); %temporary until I find the data for this
                
                % Grab IQ.
                iq(sub_count) = subs.iq(subs.subID == sub(sub_count));

                clear data_temp
                
            end % end j
            
        end % ismember
        
    end % end i
    
    save(fullfile(rootDir, ['supportFiles/devti_data_' wm{w} '_mrtrix3act.mat']), 'sub', 'age', 'sex', 'iq', 'roi', 'm')
    
    clear sub age roi m
    
end % end w



