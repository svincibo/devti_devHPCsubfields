% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

blprojectid = 'proj-5e5672430f7fa65e1d3c9621';

% Set working directories.
rootDir = '/Volumes/Seagate/devti_devHPCsubfields/';

% Read in behavioral data.
load(fullfile(rootDir, 'supportFiles/data.mat'));

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behavioral/demographic data.
    if ~isempty(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        if ~isempty(sub_contents_snr)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get standard deviation of SNR in X, Y, and Z directions.
            sd(sub_count) = std([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get training group.
            group(sub_count) = data(find((data(:, 1) == str2num(grp_contents(s).name(5:7)))), 4);
            
            clear data_snr_temp get_temp
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

% Concatenate and sort according to group.
toplot = cat(2, subID', group', m', sd', b0');
toplot = sortrows(toplot, [2 3 1], 'ascend');

% cheap trick
clear subID group m sd b0
subID = toplot(:, 1)';
group = toplot(:, 2)';
m = toplot(:, 3)';
sd = toplot(:, 4)';
b0 = toplot(:, 5)';

% SNR plot
figure(1)
hold on;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0;
alphablend = .8;

yc_color  = [50 180 100]/255; 
oc_color = [50 100 180]/255; 
a_color = [100 50 180]/255;

% yc_color = [0.6350 0.0780 0.1840]; %red
% oc_color = [0 0.4470 0.7410]; %blue
% a_color = [0.41176 0.41176 0.41176]; %gray

gscatter(m, 1:length(m), group, [yc_color; oc_color; a_color], '.', 20)
hold on;
gscatter(b0, 1:length(b0), group, [yc_color; oc_color; a_color], 'x', 8)
% plot([5 5], [0 length(subID)+0.5], ':k')

for p = 1:length(m)
    
    if group(p) == 1
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', yc_color)
        
    elseif group(p) == 2
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', oc_color)
        
    elseif group(p) == 3
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', a_color)
              
    end
    
end
gscatter(m, 1:length(m), group, [yc_color; oc_color; a_color], '.', 20)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [min([b0 m + sd])-5 max([b0 m + sd])+5];
xax.TickValues = floor(min([b0 m + sd])):10:ceil(max([b0 m + sd]));
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
% xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(subID)+0.5];
yax.TickValues = 1:length(subID);
yax.TickDirection = 'out';
yax.TickLabels = subID;
yax.FontName = fontname;
yax.FontSize = 6;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adolescents', 'Adults'}, 'Location', 'northeast');
legend box off

a.XLabel.String = 'SNR';
a.XLabel.FontSize = fontsize;
pbaspect([1 1 1])

a.YLabel.String = 'Subject ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = 'italic';

print(fullfile(rootDir, 'plots', 'plot_snr'), '-dpng')

hold off;