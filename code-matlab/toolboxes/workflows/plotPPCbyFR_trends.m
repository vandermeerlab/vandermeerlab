%% Script to plot PPC by firing rate, where trials are binned by z-scored firing-rates

cd('D:\RandomVstrAnalysis\PPCbyFR\'); % Change this to your local machine location for results
% Setting up parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
min_trial_spikes = 50; % Minimum number of spikes in a trial for it to be considered
min_trials = 25; % Minimum number of spike-count thresholded trials in a cell for it to be considered
edges = [-100, -2, -1, 0, 1, 2, 100]; % The first and last bins are arbitrarily large to contain anything outside 2 standard deviations

% Collect all FSIs
paths = FindFiles('*FSI_ppcbyfr.mat');
fsi_out = cell(length(paths), length(f_list));
for iC = 1:length(paths)
    load(paths{iC});
    for iF = 1:length(f_list)
        f_idx = find(round(od.freqs) >= f_list{iF}(1) & round(od.freqs) <= f_list{iF}(2));        
        binned_mean = zeros(1,length(edges)-1);
        for iB = 1:length(binned_mean)
            if od.binned_spk_count(iB) < min_trial_spikes
                binned_mean(iB) = NaN;
            else
                binned_mean(iB) = mean(od.binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        binned_mean(~isfinite(binned_mean)) = NaN;
        % normalize it
        isnum = ~isnan(binned_mean);
        binned_mean(isnum) = (binned_mean(isnum) - min(binned_mean(isnum))) /...
            (max(binned_mean(isnum) - min(binned_mean(isnum))));
        fsi_out{iC,iF} = binned_mean;
    end
end

% Collect all MSNs
paths = FindFiles('*MSN_ppcbyfr.mat');
msn_out = cell(length(paths), length(f_list));
for iC = 1:length(paths)
    load(paths{iC});
    for iF = 1:length(f_list)
        f_idx = find(round(od.freqs) >= f_list{iF}(1) & round(od.freqs) <= f_list{iF}(2));        
        binned_mean = zeros(1,length(edges)-1);
        for iB = 1:length(binned_mean)
            if od.binned_spk_count(iB) < min_trial_spikes
                binned_mean(iB) = NaN;
            else
                binned_mean(iB) = mean(od.binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        binned_mean(~isfinite(binned_mean)) = NaN;
        % normalize it
        isnum = ~isnan(binned_mean);
        binned_mean(isnum) = (binned_mean(isnum) - min(binned_mean(isnum))) /...
            (max(binned_mean(isnum) - min(binned_mean(isnum))));
        msn_out{iC,iF} = binned_mean;
    end
end

%% Plot trend
fig = figure('WindowState', 'maximized');
for iF = 1:length(f_list)
    % Plot FSIs
    subplot(2,length(f_list),iF)
    hold on;
    for iC = 1:length(fsi_out)
        plot(fsi_out{iC}(1,:), 'color',c_list{iF}, 'LineWidth', 0.25);
    end
    temp_ppc = cell2mat(fsi_out(:,iF));
    plot(mean(temp_ppc,'omitnan'), 'color', 'black', 'LineWidth', 3);
    xticks([1.5, 2.5, 3.5, 4.5, 5.5])
    xticklabels({'-2', '-1', '0', '1', '2'})
    xlim([0.85 6.15])
    ylabel('Norm PPC')
    xlabel('Z-scored FR')
    title(sprintf("FSI %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
    
    % Plot MSNs
    subplot(2,length(f_list),iF+length(f_list))
     hold on;
    for iC = 1:length(msn_out)
        plot(msn_out{iC}(1,:), 'color',c_list{iF}, 'LineWidth', 0.25);
    end
    temp_ppc = cell2mat(msn_out(:,iF));
    plot(mean(temp_ppc,'omitnan'), 'color', 'black', 'LineWidth', 3);
    xticks([1.5, 2.5, 3.5, 4.5, 5.5])
    xticklabels({'-2', '-1', '0', '1', '2'})
    xlim([0.85 6.15])
    ylabel('Norm PPC')
    xlabel('Z-scored FR')
    title(sprintf("MSN %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
end
%% Plot heatmap
fig = figure('WindowState', 'maximized');
for iF = 1:length(f_list)
    % Plot FSIs
    subplot(2,length(f_list),iF)
    hold on;
    temp_ppc = cell2mat(fsi_out(:,iF));
    % Sort according to second ppc bin
    [~,sidx] = sort(temp_ppc(:,2));
    temp_ppc(:,:) = temp_ppc(sidx,:);
    imagesc(temp_ppc)
    xticks([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
    xticklabels({'-2', '-1.2', '-0.4', '0.4', '1.2', '2'})
    ylabel('Norm PPC')
    xlabel('Z-scored FR')
    axis off
    title(sprintf("FSI:%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
    
    % Plot MSNs
    subplot(2,length(f_list),iF+length(f_list))
    hold on;
    temp_ppc = cell2mat(msn_out(:,iF));
    % Sort according to the second ppc bin
    [~,sidx] = sort(temp_ppc(:,2));
    temp_ppc(:,:) = temp_ppc(sidx,:);
    imagesc(temp_ppc)
    xticks([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
    xticklabels({'-2', '-1.2', '-0.4', '0.4', '1.2', '2'})
    ylabel('Norm PPC')
    xlabel('Z-scored FR')
    axis off
    title(sprintf("MSN %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
end