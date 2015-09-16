function [ pos_x,pos_y,Timestamps ] = xy_targets_shortcut( unique_folder,num_targets )
unique_id = unique_folder(1:15);

% defaults
chunksize = 5000;

% Load raw data
[Timestamps, X, Y, Angles, Targets, Points, Header] = ...
    Nlx2MatVT([unique_id,'-VT1.nvt'], [1 1 1 1 1 1], 1, 1, []);

% % Index into the data for troubleshooting
% Targets = Targets(:,1000:20000);
% Timestamps = Timestamps(:,1000:20000);

% Convert timestamps to seconds
Timestamps = Timestamps * 10^-6;

nBins = length(Targets);
chunks = [1:chunksize - 1:nBins, nBins];
pos_x = zeros(num_targets, nBins);
pos_y = zeros(num_targets, nBins);

for i = 1:length(chunks)-1
    disp(sprintf('Processing portion %d of %d', i, length(chunks)-1));
    small_target_matrix = Targets(:,chunks(i):chunks(i+1));
    [x_temp, y_temp, color_temp, valid_targets_temp] = ...
        ExtractFromTargets(small_target_matrix);
    no_data = find(x_temp == 0 & y_temp == 0);
    x_temp(no_data) = nan;
    y_temp(no_data) = nan;
    
    x_temp = x_temp(:,1:num_targets);
    y_temp = y_temp(:,1:num_targets);
    pos_x(1:num_targets, chunks(i):chunks(i)+size(x_temp,1)-1) = x_temp';
    pos_y(1:num_targets, chunks(i):chunks(i)+size(y_temp,1)-1) = y_temp';  
end
save(['pos-xytime-',unique_folder(1:15)],'pos_x','pos_y','Timestamps');
end