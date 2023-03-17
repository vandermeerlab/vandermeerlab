%% Script to plot PPC by firing rate for each cell separately

cd('D:\RandomVstrAnalysis\PPCbyFR\'); % Change this to your local machine location for results
% Setting up parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
min_trial_spikes = 50; % Minimum number of spikes in a trial for it to be considered
min_trials = 25; % Minimum number of spike-count thresholded trials in a cell for it to be considered
edges = [-100, -2, -1, 0, 1, 2, 100]; % The first and last bins are arbitrarily large to contain anything outside 2 SDs

% Plot FSIs first
paths = FindFiles('*FSI_ppcbyfr.mat');
fsi_out = cell(length(paths), length(f_list));
for iC = 1:length(paths)
    load(paths{iC});
    fig = figure('WindowState', 'maximized');
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
        
        % Plot raw PPC vs raw fr
        this_mean = od.binned_mean;
        this_sd = od.binned_sd;
        ax = subplot(3,length(f_list),iF);
        plot(binned_mean, 'color', c_list{iF}, 'LineWidth', 3);
        ax.YAxis.Exponent = 0;
        xticks([1.5, 2.5, 3.5, 4.5, 5.5])
        xticklabels({num2str(this_mean - 2*this_sd, '%.1f'), num2str(this_mean - 1*this_sd, '%.1f'), num2str(this_mean, '%.1f'), ...
            num2str(this_mean + 1*this_sd, '%.1f'), num2str(this_mean + 2*this_sd, '%.1f')});
        xlabel('Firing rate (Hz)')
        ylabel('PPC')
        title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
        
        % normalize it
        isnum = ~isnan(binned_mean);
        binned_mean(isnum) = (binned_mean(isnum) - min(binned_mean(isnum))) /...
            (max(binned_mean(isnum) - min(binned_mean(isnum))));
        
        % Plot norm PPC vs z-scored FR
        subplot(3,length(f_list),iF+length(f_list))
        plot(binned_mean, 'color', c_list{iF}, 'LineWidth', 3);
        xticks([1.5, 2.5, 3.5, 4.5, 5.5])
        xticklabels({'-2', '-1', '0', '1', '2'})
        xlabel('Z-Scored Firing rate')
        ylabel('Norm PPC')

    end

    % Load the results file to get the msn spike_count distribution
    label = strsplit(paths{iC},'_FSI');
    label = strsplit(label{1}, '\');
    label = label{end};
    toks = strsplit(label, '-');
    full_res = load(strcat('D:\RandomVstrAnalysis\final_results\', strjoin(toks(1:4),'-'),'_ft_spec.mat'));
    this_fsi = find(full_res.od.cell_type == 2);
    this_fsi = find(startsWith(full_res.od.label(this_fsi), label));

    this_freqs = od.freqs;
    this_sts = full_res.od.fsi_res.near_hfr_spec{this_fsi}.subsampled_sts;
    this_ppc = full_res.od.fsi_res.near_hfr_spec{this_fsi}.subsampled_ppc;
    
    % Plot PPC and STS for comparison
    subplot(3,length(f_list), 2*length(f_list)+1:2*length(f_list)+2)
    plot(this_freqs, this_ppc, 'color', 'black')
    hold on;
    for iF = 1:length(f_list)
        f_idx = find(round(this_freqs) >= f_list{iF}(1) & ...
            round(this_freqs) <= f_list{iF}(2));
        area(this_freqs(f_idx), this_ppc(f_idx), ...
            'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
    end
    xlabel('Frequencies (Hz)');
    ylabel('PPC');
    title('ALL trials')

    subplot(3,length(f_list), 2*length(f_list)+3:2*length(f_list)+4)
    plot(this_freqs, this_sts, 'color', 'black')
    hold on;
    for iF = 1:length(f_list)
        f_idx = find(round(this_freqs) >= f_list{iF}(1) & ...
            round(this_freqs) <= f_list{iF}(2));
        area(this_freqs(f_idx), this_sts(f_idx), ...
            'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
    end
    xlabel('Frequencies (Hz)');
    ylabel('STS');
    title('ALL trials')

    % Plot spike distribution in bins
    subplot(3,length(f_list), 3*length(f_list))
    bar(od.binned_spk_count)
    ylabel("Spike count per bin");

    sgtitle(sprintf("FSI: %s", label), 'Interpreter', 'None');
    WriteFig(fig,cat(2,'FSI_', label),1);
    close;
end   
%% Plot MSNs next
paths = FindFiles('*MSN_ppcbyfr.mat');
msn_out = cell(length(paths), length(f_list));
for iC = 1:length(paths)
    load(paths{iC});
    fig = figure('WindowState', 'maximized');
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
        
        % Plot raw PPC vs raw fr
        this_mean = od.binned_mean;
        this_sd = od.binned_sd;
        ax = subplot(3,length(f_list),iF);
        plot(binned_mean, 'color', c_list{iF}, 'LineWidth', 3);
        ax.YAxis.Exponent = 0;
        xticks([1.5, 2.5, 3.5, 4.5, 5.5])
        xticklabels({num2str(this_mean - 2*this_sd, '%.1f'), num2str(this_mean - 1*this_sd, '%.1f'), num2str(this_mean, '%.1f'), ...
            num2str(this_mean + 1*this_sd, '%.1f'), num2str(this_mean + 2*this_sd, '%.1f')});
        xlabel('Firing rate (Hz)')
        ylabel('PPC')
        title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
        
        % normalize it
        isnum = ~isnan(binned_mean);
        binned_mean(isnum) = (binned_mean(isnum) - min(binned_mean(isnum))) /...
            (max(binned_mean(isnum) - min(binned_mean(isnum))));
        
        % Plot norm PPC vs z-scored FR
        subplot(3,length(f_list),iF+length(f_list))
        plot(binned_mean, 'color', c_list{iF}, 'LineWidth', 3);
        xticks([1.5, 2.5, 3.5, 4.5, 5.5])
        xticklabels({'-2', '-1', '0', '1', '2'})
        xlabel('Z-Scored Firing rate')
        ylabel('Norm PPC')

    end

    % Load the results file to get the msn spike_count distribution
    label = strsplit(paths{iC},'_MSN');
    label = strsplit(label{1}, '\');
    label = label{end};
    toks = strsplit(label, '-');
    full_res = load(strcat('D:\RandomVstrAnalysis\final_results\', strjoin(toks(1:4),'-'),'_ft_spec.mat'));
    this_msn = find(full_res.od.cell_type == 1);
    this_msn = find(startsWith(full_res.od.label(this_msn), label));

    this_freqs = od.freqs;
    this_sts = full_res.od.msn_res.near_hfr_spec{this_msn}.sts_vals;
    this_ppc = full_res.od.msn_res.near_hfr_spec{this_msn}.ppc;
    
    % Plot PPC and STS for comparison
    subplot(3,length(f_list), 2*length(f_list)+1:2*length(f_list)+2)
    plot(this_freqs, this_ppc, 'color', 'black')
    hold on;
    for iF = 1:length(f_list)
        f_idx = find(round(this_freqs) >= f_list{iF}(1) & ...
            round(this_freqs) <= f_list{iF}(2));
        area(this_freqs(f_idx), this_ppc(f_idx), ...
            'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
    end
    xlabel('Frequencies (Hz)');
    ylabel('PPC');
    title('ALL trials')

    subplot(3,length(f_list), 2*length(f_list)+3:2*length(f_list)+4)
    plot(this_freqs, this_sts, 'color', 'black')
    hold on;
    for iF = 1:length(f_list)
        f_idx = find(round(this_freqs) >= f_list{iF}(1) & ...
            round(this_freqs) <= f_list{iF}(2));
        area(this_freqs(f_idx), this_sts(f_idx), ...
            'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
    end
    xlabel('Frequencies (Hz)');
    ylabel('STS');
    title('ALL trials')

    % Plot spike distribution in bins
    subplot(3,length(f_list), 3*length(f_list))
    bar(od.binned_spk_count)
    ylabel("Spike count per bin");

    sgtitle(sprintf("MSN: %s", label), 'Interpreter', 'None');
    WriteFig(fig,cat(2,'MSN_', label),1);
    close;
end