function [ x,y,Timestamps ] = xy_targets_chunks( unique_folder,chunksize,var_tolerance )
unique_id = unique_folder(1:15);

[Timestamps, X, Y, Angles, Targets, Points, Header] = ...
    Nlx2MatVT([unique_id,'-VT1.nvt'], [1 1 1 1 1 1], 1, 1, []);

% Index into the data for troubleshooting
% Targets = Targets(:,5000:45000);
% Timestamps = Timestamps(:,5000:45000);

nBins = length(Targets);
chunks = [1:chunksize-1:nBins, nBins];
x = zeros(nBins,1);
y = zeros(nBins,1);
current = 1;
total = ceil(nBins/(chunksize-1));

for i = 1:length(chunks)-1
    disp(sprintf('Processing portion %d of %d',current,total));
    small_target_matrix = Targets(:,chunks(i):chunks(i+1));
    [x_temp, y_temp, color_temp, valid_targets_temp] = ...
        ExtractFromTargets(small_target_matrix);
    x_temp(x_temp == 0) = nan;
    y_temp(y_temp == 0) = nan;
        
    % Averaging over targets when varience isn't too high, otherwise nan.
    x_avg = nanmean(x_temp');
    x_var = nanvar(x_temp');
    x_avg(x_var > var_tolerance) = nan;
    x(chunks(i):chunks(i)+length(x_temp(:,1))-1,:) = x_avg;
    
    y_avg = nanmean(y_temp');
    y_var = nanvar(y_temp');
    y_avg(y_var > var_tolerance) = nan;
    y(chunks(i):chunks(i)+length(y_temp(:,1))-1,:) = y_avg;
    
%     % Using the first target only.
%     x(chunks(i):chunks(i)+length(x_temp(:,1))-1,:) = x_temp(:,1);
%     y(chunks(i):chunks(i)+length(y_temp(:,1))-1,:) = y_temp(:,1);
    
    current = current + 1;
end
end

