% This script reads in the tract profiles from Brad Caron's Tract Profile
% App and makes them into a group-level tract-profile.

clear all; close all; clc
format long g

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Select WM measure.
wm = {'fa', 'md'};

% Set working directories.
rootdir = '/Volumes/Seagate/devti_devHPCsubfields';

% Get age.
sex_in = load(fullfile(rootdir, 'supportFiles', 'devti_sex_all.mat'));
age_in = load(fullfile(rootdir, 'supportFiles', 'age_n=79.mat'));
task1 = load(fullfile(rootdir, 'taskperf', 'task1_behavior_n=78.mat'));
t1 = array2table(cat(2, task1.subnumbers', task1.allPerfData), 'VariableNames', ['subnumbers' task1.allPerfDataCols]);
% Keep only columsn of interest.
beh = t1(:, [1 13:16]); 
% Rename to help my brain.
beh.Properties.VariableNames{1} = 'subID';
beh.Properties.VariableNames{2} = 'assoc';
beh.Properties.VariableNames{3} = 'infer';
beh.Properties.VariableNames{4} = 'infassoc';
% t3 = array2table(task3.allPerfData, 'VariableNames', task3.allPerfDataCols);
% task4 = load(fullfile(rootdir, 'taskperf', 'task4_behavior_n=78.mat'));
% task5 = load(fullfile(rootdir, 'taskperf', 'task5_behavior_n=78.mat'));

%% ROI.

% Get contents of the directory where the subject MRI data are stored.
grp_contents = dir(fullfile(rootdir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Do for each wm measure.
for w = 1:size(wm, 2)
    
    if strcmp(wm{w}, 'fa')
        % Select roi data to collect.
        rois = {'fa_dwi_ROI_b_hip.nii.gz', 'fa_dwi_ROI_b_body.nii.gz', 'fa_dwi_ROI_b_head.nii.gz', 'fa_dwi_ROI_b_tail.nii.gz', ...
            'fa_dwi_ROI_b_ca1.nii.gz', 'fa_dwi_ROI_b_ca23.nii.gz', 'fa_dwi_ROI_b_dg.nii.gz', 'fa_dwi_ROI_b_sub.nii.gz', ...
            'ca1_head', 'ca23_head', 'dg_head', 'sub_head', ...
            'ca1_body', 'ca23_body', 'dg_body', 'sub_body'};
    elseif strcmp(wm{w}, 'md')
        % Select roi data to collect.
        rois = {'md_dwi_ROI_b_hip.nii.gz', 'md_dwi_ROI_b_body.nii.gz', 'md_dwi_ROI_b_head.nii.gz', 'md_dwi_ROI_b_tail.nii.gz', ...
            'md_dwi_ROI_b_ca1.nii.gz', 'md_dwi_ROI_b_ca23.nii.gz', 'md_dwi_ROI_b_dg.nii.gz', 'md_dwi_ROI_b_sub.nii.gz', ...
            'ca1_head', 'ca23_head', 'dg_head', 'sub_head', ...
            'ca1_body', 'ca23_body', 'dg_body', 'sub_body'};
    end
    
    % Load in each ROI's measures for this subject.
    sub_count = 0;
    for s = 1:size(grp_contents, 1)
        
        % Only read in data for subjects that are in the age file and also have diffusion data.
        if sum(ismember(age_in.subnumbers, grp_contents(s).name(end-2:end))) ~= 0
            
            disp(grp_contents(s).name(end-2:end))
            
            sub_count = sub_count +1;
            
            % Grab subID.
            sub_temp = str2num(grp_contents(s).name(end-2:end));
            
            % Grab age.
            age_cov_temp = age_in.day2age(find(age_in.subnumbers == sub_temp));
            
            % Assign age group.
            if age_cov_temp >= 18.0
                
                age_group_temp = 3; % adult
                
            elseif age_cov_temp < 12.0
                
                age_group_temp = 1; % child
                
            else
                
                age_group_temp = 2; % adolescent
                
            end
            
            % Grab sex. 1 = M, 2 = F
            sex_temp = sex_in.sex(find(sex_in.subID == sub_temp));
            
            % % Grab performance.
            assoc = beh.assoc(beh.subID == sub_temp);
            infer = beh.infer(beh.subID == sub_temp);
         
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_rois = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_*.nii.gz']));
            
            % Remove the '.' and '..' files.
            sub_contents_rois = sub_contents_rois(arrayfun(@(x) x.name(1), sub_contents_rois) ~= '.');
            
            % cycle through rois
            for r = 1:size(rois, 2)
                
                % Find where this roi is in the folder.
                idx = find(contains({sub_contents_rois.name}, rois{r}));
                
                if ~isempty(idx) && r < 9
                    
                    % Read in data for this subject and this roi.
                    data_temp = niftiread(fullfile(sub_contents_rois(idx).folder, sub_contents_rois(idx).name));
                    
                    % Convert all zeros to NaN.
                    data_temp(data_temp == 0) = NaN;
                    
                    % Average and standard deviation all non-NaN values.
                    m_temp = nanmean(nanmean(nanmean(data_temp, 1), 2), 3);
                    %sd{sub_count, j} = nanstd(nanstd(nanstd(data_temp, 1), 2), 3);
                    
                    % Grab roi name.
                    roi_temp = sub_contents_rois(idx).name(14:end-7);
                    
                elseif ~isempty(idx)
                    
                    
                    % This is one of the conjunction ROIs.
                    if strcmp(rois{r}, 'ca1_head')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_ca1.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_head.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'ca23_head')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_ca23.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_head.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'dg_head')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_dg.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_head.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'sub_head')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_sub.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_head.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'ca1_body')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_ca1.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_body.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'ca23_body')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_ca23.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_body.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'dg_body')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_dg.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_body.nii.gz']));
                        second = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                    elseif strcmp(rois{r}, 'sub_body')
                        
                        % Read in first roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_sub.nii.gz']));
                        first = niftiread(fullfile(conjroi(1).folder, conjroi(1).name));
                        clear conjroi;
                        
                        % Read in second roi.
                        conjroi = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, ['/dt-raw.tag-tensor_metrics*/' wm{w} '*_ROI_b_body.nii.gz']));
                        
                    end
                    
                    % Find the union of these two ROIs.
                    idx1 = first ~= 0; idx2 = second ~= 0;
                    new_roi = idx1.*idx2;
                    
                    % Get new_roi with measures.
                    data_temp = new_roi.*first;
                    
                    % Convert all zeros to NaN.
                    data_temp(data_temp == 0) = NaN;
                    
                    % Return the mean microstructure measure for the conjunction roi.
                    m_temp = nanmean(nanmean(nanmean(data_temp, 1), 2), 3);
                    
                    % Grab roi name.
                    roi_temp = rois{r};
                    
                end % ifexist
                
                if s == 1 && r == 1
                    
                    m = table(sub_temp, age_cov_temp, age_group_temp, sex_temp, assoc, infer, m_temp, 'VariableNames', {'subID', 'age_cov', 'age_group', 'sex', 'assoc', 'infer', roi_temp});
                    
                elseif s == 1
                    
                    m = cat(2, m, table(m_temp));
                    m.Properties.VariableNames{end} = roi_temp;
                    
                end
                
                if s ~= 1
                    
                    if r == 1
                        % append demo values
                        m(end+1, 1:6) = table(sub_temp, age_cov_temp, age_group_temp, sex_temp, assoc, infer);
                    end
                    
                    % find row idx
                    s_idx = find(m.subID == sub_temp);
                    
                    % find column idx
                    r_idx = find(strcmp(m.Properties.VariableNames, roi_temp));
                    
                    m(s_idx, r_idx) = table(m_temp);
                    
                end
                
                clear data_temp
                
            end % end r
            
        end % ismember
        
    end % end s
    
    filename = sprintf('devti_data_%s_%s', wm{w}, datestr(now,'yyyymmdd'));
    
    % save it as a matlab table
    save(fullfile(rootdir, 'supportFiles', filename), 'm')
    
    % Save as a CSV files.
    writetable(m, fullfile(rootdir, 'supportFiles', [filename '.csv']))
    
    clear m
    
end % end w



