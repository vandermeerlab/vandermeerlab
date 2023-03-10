%% Script to plot trialwise PPC summary

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
msn_sd = [];
fsi_sd = [];
min_trial_spikes = 50; % Minimum number of spikes per trial
min_trials = 2; % Minimum number of trials in a cell
min_freq = 2; % Minimum frequency in Hertz

for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        % do fsi stuff
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ~od.fsi_res.near_spec{iC}.flag_no_control_split
                % do fsi_stuff
                nz_trials = find(od.fsi_res.near_spec{iC}.trialwise_spk_count >= min_trial_spikes);
                if length(nz_trials) < min_trials
                    continue;
                end
                nz2 = find(od.fsi_res.near_spec{iC}.trialwise_spk_count > 0);
                x1 = find(round(od.fsi_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_ppc = od.fsi_res.near_spec{iC}.trialwise_ppc(nz_trials,x1:end);
                fsi_sd = [fsi_sd; std(this_ppc)];
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                nz_trials = find(od.msn_res.near_spec{iC}.trialwise_spk_count >= min_trial_spikes);
                if length(nz_trials) < min_trials
                    continue;
                end
                nz2 = find(od.msn_res.near_spec{iC}.trialwise_spk_count > 0);
                x1 = find(round(od.msn_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_ppc = od.msn_res.near_spec{iC}.trialwise_ppc(nz_trials,x1:end);
                msn_sd = [msn_sd; std(this_ppc)];
            end
        end
    end
end
%% Plot histograms
close all;
q0 = mean(fsi_sd, 2);
q1 = mean(msn_sd, 2);
fig = figure('WindowState', 'maximized');
h1 = histogram(q0,0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(q1,0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax = gca(fig);
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
ax.TickDir = 'out';
ax.YTick = [0 0.5 1];
ax.XTick = [0 0.3 0.5];
ax.XLabel.String = 'Mean of SDs over all frequencies';
ax.YLabel.String = 'Proportion';
% ax.XLabel.FontSize = 40;
ax.XLim = [0, 0.5];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 40;
leg.FontWeight = 'bold';
box off;
% title('Distribution of Average SD', 'FontSize', 40)

