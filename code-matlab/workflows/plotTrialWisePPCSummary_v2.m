%% Script to plot trialwise PPC summary 

cd('D:\RandomVstrAnalysis\temp2\'); % Change this to your local machine location for results
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
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'num_subsamples') && od.fsi_res.near_spec{iC}.num_subsamples~=0
                x1 = find(round(od.fsi_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_ppc = od.fsi_res.near_spec{iC}.trialwise_ppc(:,:,x1:end);            
                % Method 1 (old method): The PPC for each trial is the average over all the subsamples for that trial
                tw_ppc = squeeze(nanmean(this_ppc,1)); % average ppc for each frequency in each trial, over subsamples
                sd1 = nanstd(tw_ppc,1,1); % standard deviation for each frequency, over trials
                sd1 = mean(sd1); % average sd

                % Method 2: (Find average SD for each subsample and then average it across subsamples
                sd2 = squeeze(nanstd(this_ppc,1,2)); % standard deviation for each frquency across trials, in each of the subsamples
                sd2 = mean(sd2,1); % 'average sd' for each subsample
                sd2 = mean(sd2); % 'average sd' averaged over all subsamples

                % Method 3: (Save the histograms? ask mvdm)

                fsi_sd = [fsi_sd; sd1, sd2];
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'trialwise_ppc')
                % do msn_stuff
                x1 = find(round(od.msn_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_ppc = od.msn_res.near_spec{iC}.trialwise_ppc(:,x1:end); 
                sd = nanstd(this_ppc,1,1); % standard deviation for each frequency, over trials
                sd = mean(sd); % average sd
                msn_sd = [msn_sd; sd];
            end
        end
    end
end
%% Plot histograms
fig = figure('WindowState', 'maximized');
ax1 = subplot(2,1,1);
h1 = histogram(fsi_sd(:,1), 0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(msn_sd,0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.TickDir = 'out';
ax1.YTick = [0 0.5 1];
ax1.XTick = [0:0.1:0.5];
ax1.XLabel.String = 'Mean of SDs over all frequencies';
ax1.YLabel.String = 'Proportion';
ax1.XLim = [0 0.5];
ax1.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax1.Title.String = 'Old method';
ax1.Box = 'off';


ax2 = subplot(2,1,2);
h1 = histogram(fsi_sd(:,2), 0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(msn_sd,0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax2.XAxis.FontSize = 20;
ax2.YAxis.FontSize = 20;
ax2.TickDir = 'out';
ax2.YTick = [0 0.5 1];
ax2.XTick = [0:0.1:0.5];
ax2.XLabel.String = 'Mean of SDs over all frequencies';
ax2.YLabel.String = 'Proportion';
ax2.XLim = [0, 0.5];
ax2.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax2.Title.String= 'New method';


