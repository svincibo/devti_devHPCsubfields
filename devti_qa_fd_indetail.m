clear all; close all; clc
format shortG

yc_color  = [50 180 100]/255; 
oc_color = [50 100 180]/255; 
a_color = [100 50 180]/255;

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/Seagate/devti_devHPCsubfields/';
% addpath(genpath(fullfile(rootDir, 'proj-5e5672430f7fa65e1d3c9621')));

remove_outliers = 'yes';
include = 'all'; % all, childrenonly
outlier = [11 90];

%%%%%%%%%%%%%%% TESTING AREA %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/devti_data_beh_forSPSS_20220705.csv'], 'TreatAsEmpty', {'.', 'na'});

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
            % Remove the '.' and '..' files.
            sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
            
            % Get motion parameters for this subject.
            load(fullfile(sub_contents_motion.folder,sub_contents_motion.name));
            
            % Select only the translation/rotation parameters.
            mot = vertcat(xform(:).ecParams);
            mot = mot(:, 1:6); % xyz translation xyz rotation
            
            % Motion parameters represent the translation/rotation that must occur
            % to return the image at timepoint tn to the place that it was at timepoint
            % t0. Thus, to calculate instantaneouse parameters, we need a moving
            % difference.
            for m = 1:size(mot, 2)
                
                % Get moving difference for each parameter. Append row of zeros for t1, per convention (Power et al., 2014).
                % This step creates the fd timecourses for each motion parameter for each run in the motion struct and represents instantaneous movement.
                movingdifference(:, m) = [0 ; diff(mot(:, m), 1, 1)]';
                
            end
            
            % Get an overall fd for all 6 parameters for each run.
            % This step creates the fd summary statistic for all 6 parameters for each timepoint in a run for each run (e.g., scalar FD timecourse).
            motion(sub_count, :) = sum(abs(movingdifference), 2)';
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name(5:7);
            
            % Get training group.
            group(sub_count) = beh_data_in_tbl.group(s);  
    
end % end s

meanmotion = mean(motion, 2);
sd = std(motion, 0, 2);

% Concatenate and sort according to group.
toplot = cat(2, subID', group', meanmotion, sd);
toplot = sortrows(toplot, [2 3 1], 'ascend');

% FD plot
figure(1)
hold on;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
xticklength = 0;
alphablend = .8;

gscatter(toplot(:, 3), 1:length(toplot(:, 3)), toplot(:, 2), [yc_color; oc_color; a_color], '.', 30)

hold on;
for p = 1:length(toplot(:, 3))

    if toplot(p, 2) == 1
        color = yc_color;
    elseif toplot(p, 2) == 2
        color = oc_color;
    elseif toplot(p, 2) == 3
        color = a_color;
    end
    plot([toplot(p, 3) - abs(toplot(p, 4)) toplot(p, 3) + abs(toplot(p, 4))], [p p], 'Color', color)

end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 2];
xax.TickValues = 0:0.5:2;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(lab)+0.5];
yax.TickValues = 1.5:2:ceil(length(lab));
yax.TickDirection = 'out';
yax.TickLabels = lab;
yax.FontName = fontname;
yax.FontSize = 8;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adolescents', 'Adults'}, 'Location', 'northeast');
legend box off

a.XLabel.String = 'Framewise Displacement (FD)';
a.XLabel.FontSize = fontsize;

a.YLabel.String = 'Subject ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';

pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_fd'), '-dpng')

hold off;