%% Script to plot trialwise PPC summary 

cd('D:\RandomVstrAnalysis\trialwise_ppc\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
fsi_old_sd = [];
fsi_sd = [];
msn_sd = [];
fsi_sess = {};
msn_sess = {};
min_trial_spikes = 50; % Minimum number of spikes per trial
min_trials = 25; % Minimum number of trials in a cell
min_freq = 2; % Minimum frequency in Hertz
fsi_sp = {};
msn_sp = {};
% Do this only for significantly phase_locked neurons
load('D:\RandomVstrAnalysis\final_results\significantly_phase_locked.mat');

for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        fsi_sess{jdx} = [];
        msn_sess{jdx} = [];
        fsi_sp{jdx} = [];
        msn_sp{jdx} = [];
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'num_subsamples') && od.fsi_res.near_spec{iC}.num_subsamples~=0
                for iS = 1:length(sig3_fsi)
                    sig_idx = 0;
                    this_label = sig3_fsi{iS,1};
                    this_label = this_label{1};
                    if strcmp(this_label, fsi_labels{iC})
                        sig_idx = iS;
                        break
                    end
                end             
                if (sig_idx == 0) || ~any(sig3_fsi{sig_idx,2})
                    continue;
                end
                x1 = find(round(od.fsi_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_ppc = od.fsi_res.near_spec{iC}.trialwise_ppc(:,:,x1:end);            
                % Method 1 (old method): The PPC for each trial is the average over all the subsamples for that trial
                tw_ppc = squeeze(nanmean(this_ppc,1)); % average ppc for each frequency in each trial, over subsamples
                valid_trials = ~isnan(tw_ppc(:,1));
                if sum(valid_trials) >= min_trials
                    sd1 = std(tw_ppc(valid_trials,:),1,1); % standard deviation for each frequency, over trials
                    sd1 = mean(sd1); % average sd
                    fsi_old_sd = [fsi_old_sd; sd1];

                    % Method 2: (Find average SD for each subsample and then average it across subsamples
                    sd2 = squeeze(nanstd(this_ppc,1,2)); % standard deviation for each frquency across trials, in each of the subsamples
                    sd2 = mean(sd2,2); % 'average sd' for each subsample
        %                 sd2 = mean(sd2); % 'average sd' averaged over all subsamples
                    
                    fsi_sd = [fsi_sd; sd2];
                    fsi_sess{jdx} = [fsi_sess{jdx}; sd2];
                    this_sum = sum(this_ppc(:,valid_trials,:), 3); %sum of ppc values over all frequencies, for each subsample, and each trial
                    this_sum = mean(this_sum,2);
                    fsi_sp{jdx} = [fsi_sp{jdx};this_sum];
                end
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'trialwise_ppc')
                for iS = 1:length(sig3_msn)
                    sig_idx = 0;
                    this_label = sig3_msn{iS,1};
                    this_label = this_label{1};
                    if strcmp(this_label, msn_labels{iC})
                        sig_idx = iS;
                        break
                    end
                end             
                if (sig_idx == 0) || ~any(sig3_msn{sig_idx,2})
                    continue;
                end
                % do msn_stuff
                x1 = find(round(od.msn_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_ppc = od.msn_res.near_spec{iC}.trialwise_ppc(:,x1:end);
                valid_trials = ~isnan(this_ppc(:,1));
                if sum(valid_trials) >= min_trials 
                    sd = std(this_ppc(valid_trials,:),1,1); % standard deviation for each frequency, over trials
                    sd = mean(sd); % average sd
                    msn_sd = [msn_sd; sd];
                    msn_sess{jdx} = [msn_sess{jdx}; sd];
                    this_sum = sum(this_ppc(valid_trials,:),2);
                    this_sum = mean(this_sum);
                    msn_sp{jdx} = [msn_sp{jdx}; this_sum];
                end
            end
        end
    end
end
%% Plot histograms
fig = figure('WindowState', 'maximized');
ax1 = subplot(2,1,1);
h1 = histogram(fsi_old_sd, 0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
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
h1 = histogram(fsi_sd, 0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
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
%% Look at session-wise distributions
for iS = 1:length(fsi_sess)
    if ~isempty(fsi_sess{iS})
        figure('WindowState', 'maximized');
        subplot(2,1,1)
        h1 = histogram(fsi_sess{iS}, 0:0.001:0.1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
        hold on;
        h2 = histogram(msn_sess{iS},0:0.001:0.1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
        xlim([0 0.1])
        subplot(2,1,2)
        h1 = histogram(fsi_sp{iS},-0.5:0.1:0.5, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
        hold on;
        h2 = histogram(msn_sp{iS},-0.5:0.1:0.5, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
        dummy = 1;
    end
        
end

%%
figure('WindowState', 'maximized');
for iS = 1:length(fsi_sess)
    if ~isempty(fsi_sess{iS})
        hold on;
        scatter(fsi_sp{iS},fsi_sess{iS}, 'green');
        scatter(msn_sp{iS},msn_sess{iS}, 'red');
        xlabel('Sum of PPC of all frequencies')
        ylabel('Average SD')
        dummy = 1;
    end
        
end
