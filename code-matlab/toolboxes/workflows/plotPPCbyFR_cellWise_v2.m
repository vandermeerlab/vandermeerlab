%% Script to plot PPC by firing rate for each cell separately

cd('D:\RandomVstrAnalysis\PPCbyFR\'); % Change this to your local machine location for results
% Setting up parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
min_trial_spikes = 1; % Minimum number of spikes in a trial for it to be considered
min_trials = 25; % Minimum number of spike-count thresholded trials in a cell for it to be considered
z_threshold = 3;

%% Plot MSNs first
paths = FindFiles('*MSN_ppcbyfr.mat');
msn_out = cell(length(paths), length(f_list));
for iC = 1:length(paths)
    load(paths{iC});
    if any(od.z_binned_spk_count(2:5) < min_trial_spikes)
        continue
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
    this_ppc = full_res.od.msn_res.near_hfr_spec{this_msn}.ppc';
    
    % This block will probably have to be changed after recalculating
    % values, for now let's assume norm distribution of values
    shuf_mean = mean(od.shuf_ppc);
    shuf_sd = std(od.shuf_ppc);
    z_ppc = (this_ppc - shuf_mean)./(shuf_sd);
    sig_ppc = z_ppc > z_threshold;

    % Make a list of significant frequency values in the frequency range of
    % interest
    sig_freq_idx = cell(size(f_list));  
    for iF = 1:length(f_list)
        f_idx = find(round(od.freqs) >= f_list{iF}(1) & round(od.freqs) <= f_list{iF}(2)); 
        if any(sig_ppc(f_idx))
            sig_freq_idx{iF} = f_idx(sig_ppc(f_idx));
        end    
    end
    sig_band = find(cellfun(@length, sig_freq_idx));
    
    % No need to plot if no significant values
    if isempty(sig_band) 
        continue;
    end

    fig = figure('WindowState', 'maximized');
    for iF = 1:length(sig_band) 
        
        % Plot Z-norm bin on top
        f_idx = sig_freq_idx{sig_band(iF)};
        z_mean = zeros(1, length(od.z_edges)-1);
        for iB = 1:length(z_mean)
            if od.z_binned_spk_count(iB) < min_trial_spikes
                z_mean(iB) = NaN;
            else
                z_mean(iB) = mean(od.z_binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        z_mean(~isfinite(z_mean)) = NaN;
        this_mean = od.z_binned_mean;
        this_sd = od.z_binned_sd;
        ax = subplot(3,length(f_list),sig_band(iF));
        plot(z_mean, 'color', c_list{sig_band(iF)}, 'LineWidth', 1.5);
        ax.YAxis.Exponent = 0;
        xticks([1.5, 2.5, 3.5, 4.5, 5.5])
        xticklabels({num2str(this_mean - 2*this_sd, '%.1f'),  num2str(this_mean - 1*this_sd, '%.1f'), ...
            num2str(this_mean, '%.1f'), num2str(this_mean + 1*this_sd, '%.1f'),  ...
            num2str(this_mean + 2*this_sd, '%.1f')});
        xlim([1 6])
        xlabel('Firing rate (Hz)')
        ylabel('PPC')
        title(sprintf("Z-norm:%d Hz - %d Hz", f_list{sig_band(iF)}(1), f_list{sig_band(iF)}(2)), 'color', c_list{sig_band(iF)})

        % Plot 0-1 norm next
        f_idx = sig_freq_idx{sig_band(iF)};
        n_mean = zeros(1, length(od.n_edges)-1);
        for iB = 1:length(n_mean)
            if od.n_binned_spk_count(iB) < min_trial_spikes
                n_mean(iB) = NaN;
            else
                n_mean(iB) = mean(od.n_binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        n_mean(~isfinite(z_mean)) = NaN;
        ax = subplot(3,length(f_list),sig_band(iF)+length(f_list));
        plot(n_mean, 'color', c_list{sig_band(iF)}, 'LineWidth', 1.5);
        ax.YAxis.Exponent = 0;
        xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
        % Round-about way to extract
        x_labels = {};
        for iL = 1:length(od.n_edges)-1
            x_min = min(od.mfr(od.mfr>0)); 
            x_max = max(od.mfr(od.mfr>0));
            x_labels{iL} = num2str(x_min + mean(od.n_edges(iL:iL+1))*(x_max - x_min), '%.1f');
        end
        xlim([1 6])
        xticklabels(x_labels);
        xlabel('Firing rate (Hz)')
        ylabel('PPC')
        title(sprintf("0-1 Norm:%d Hz - %d Hz", f_list{sig_band(iF)}(1), f_list{sig_band(iF)}(2)), 'color', c_list{sig_band(iF)})
    end
    
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
    ax = subplot(3,length(f_list), 2*length(f_list)+5);
    ax.YAxis.Exponent = 0;
    hold on;
    bar(1:6,[od.z_binned_spk_count';od.n_binned_spk_count'], 'BarWidth', 1.5)
    yline(50, '-black');
    ylabel("Spike count per bin");
    xticks(1:6)
    legend({'Z', '0-1'},'FontSize', 10, 'Location', 'northwest','Box','off')
    sgtitle(sprintf("MSN: %s", label), 'Interpreter', 'None');
    WriteFig(fig,cat(2,'MSN_', label),1);
    close;
end

%% Plot FSIs next
paths = FindFiles('*FSI_ppcbyfr.mat');
fsi_out = cell(length(paths), length(f_list));
for iC = 1:length(paths)
    load(paths{iC});
    if any(od.z_binned_spk_count(2:5) < min_trial_spikes)
        continue
    end

    % Load the results file to get the fsi spike_count distribution
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
    
    % This block will probably have to be changed after recalculating
    % values, for now let's assume norm distribution of values
    shuf_mean = mean(od.shuf_ppc);
    shuf_sd = std(od.shuf_ppc);
    z_ppc = (this_ppc - shuf_mean)./(shuf_sd);
    sig_ppc = z_ppc > z_threshold;

    % Make a list of significant frequency values in the frequency range of
    % interest
    sig_freq_idx = cell(size(f_list));  
    for iF = 1:length(f_list)
        f_idx = find(round(od.freqs) >= f_list{iF}(1) & round(od.freqs) <= f_list{iF}(2)); 
        if any(sig_ppc(f_idx))
            sig_freq_idx{iF} = f_idx(sig_ppc(f_idx));
        end    
    end
    sig_band = find(cellfun(@length, sig_freq_idx));
    
    % No need to plot if no significant values
    if isempty(sig_band) 
        continue;
    end

    fig = figure('WindowState', 'maximized');
    for iF = 1:length(sig_band) 
        
        % Plot Z-norm bin on top
        f_idx = sig_freq_idx{sig_band(iF)};
        z_mean = zeros(1, length(od.z_edges)-1);
        for iB = 1:length(z_mean)
            if od.z_binned_spk_count(iB) < min_trial_spikes
                z_mean(iB) = NaN;
            else
                z_mean(iB) = mean(od.z_binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        z_mean(~isfinite(z_mean)) = NaN;
        this_mean = od.z_binned_mean;
        this_sd = od.z_binned_sd;
        ax = subplot(3,length(f_list),sig_band(iF));
        plot(z_mean, 'color', c_list{sig_band(iF)}, 'LineWidth', 1.5);
        ax.YAxis.Exponent = 0;
        xticks([1.5, 2.5, 3.5, 4.5, 5.5])
        xticklabels({num2str(this_mean - 2*this_sd, '%.1f'),  num2str(this_mean - 1*this_sd, '%.1f'), ...
            num2str(this_mean, '%.1f'), num2str(this_mean + 1*this_sd, '%.1f'),  ...
            num2str(this_mean + 2*this_sd, '%.1f')});
        xlim([1 6])
        xlabel('Firing rate (Hz)')
        ylabel('PPC')
        title(sprintf("Z-norm:%d Hz - %d Hz", f_list{sig_band(iF)}(1), f_list{sig_band(iF)}(2)), 'color', c_list{sig_band(iF)})

        % Plot 0-1 norm next
        f_idx = sig_freq_idx{sig_band(iF)};
        n_mean = zeros(1, length(od.n_edges)-1);
        for iB = 1:length(n_mean)
            if od.n_binned_spk_count(iB) < min_trial_spikes
                n_mean(iB) = NaN;
            else
                n_mean(iB) = mean(od.n_binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        n_mean(~isfinite(z_mean)) = NaN;
        ax = subplot(3,length(f_list),sig_band(iF)+length(f_list));
        plot(n_mean, 'color', c_list{sig_band(iF)}, 'LineWidth', 1.5);
        ax.YAxis.Exponent = 0;
        xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
        % Round-about way to extract
        x_labels = {};
        for iL = 1:length(od.n_edges)-1
            x_min = min(od.mfr(od.mfr>0)); 
            x_max = max(od.mfr(od.mfr>0));
            x_labels{iL} = num2str(x_min + mean(od.n_edges(iL:iL+1))*(x_max - x_min), '%.1f');
        end
        xlim([1 6])
        xticklabels(x_labels);
        xlabel('Firing rate (Hz)')
        ylabel('PPC')
        title(sprintf("0-1 Norm:%d Hz - %d Hz", f_list{sig_band(iF)}(1), f_list{sig_band(iF)}(2)), 'color', c_list{sig_band(iF)})
    end
    
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
    ax = subplot(3,length(f_list), 2*length(f_list)+5);
    ax.YAxis.Exponent = 0;
    hold on;
    bar(1:6,[od.z_binned_spk_count';od.n_binned_spk_count'], 'BarWidth', 1.5)
    yline(50, '-black');
    ylabel("Spike count per bin");
    xticks(1:6)
    legend({'Z', '0-1'},'FontSize', 10, 'Location', 'northwest','Box','off')
    sgtitle(sprintf("FSI: %s", label), 'Interpreter', 'None');
    WriteFig(fig,cat(2,'FSI_', label),1);
    close;
end