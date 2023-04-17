% Plotting the distribution of variance in trial-wise normalized firing-rates

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
fsi_frv = []; msn_frv = [];                 % no thresholding, only non-nan ppc trials
fsi_frv_nt = []; msn_frv_nt = [];           % thresholding, only non-nan ppc trials
sig_fsi_frv = []; sig_msn_frv = [];         % cells with significant phase locking in at least one of the frequency bands, non thresholding, non-nan ppc trials
sig_fsi_frv_nt = []; sig_msn_frv_nt = [];   % cells with significant phase locking in at least one of the frequency bands, thresholding, non-nan ppc trials
min_trial_spikes = 50;
min_trials = 25;
% Loading the significance list
load('D:\RandomVstrAnalysis\final_results\significantly_phase_locked.mat');
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
                clean_fsi = clean_fsi + 1;
                tw_mfr = od.fsi_res.near_spec{iC}.mfr;
                sig_idx = strcmp(sig3_fsi(:,1), fsi_labels{iC});
                % Including trials which have min_trial_spikes and have no NaNs
                valid_nothresh = ~isnan(od.fsi_res.near_spec{iC}.no_thresh_trialwise_ppc(1,:,1))'; 
                nonan_thresh = ~isnan(od.fsi_res.near_spec{iC}.thresh_trialwise_ppc(1,:,1))';
                spike_thresh = od.fsi_res.near_spec{iC}.thresh_trialwise_spk_count >= min_trial_spikes;
                valid_thresh = nonan_thresh & spike_thresh;
               
                if sum(valid_thresh) >= min_trials
                    mfr_t = tw_mfr(valid_thresh);
                    mfr_t = (mfr_t - min(mfr_t))/(max(mfr_t) - min(mfr_t)); % Normalize values between 0-1
                    fsi_frv = [fsi_frv, std(mfr_t)];
                    if any(sig3_fsi{sig_idx,2})
                        sig_fsi_frv = [sig_fsi_frv, std(mfr_t)];
                    end
                end
                
                mfr_nt= tw_mfr(valid_nothresh);
                mfr_nt = (mfr_nt - min(mfr_nt))/(max(mfr_nt) - min(mfr_nt)); % Normalize values between 0-1
                fsi_frv_nt = [fsi_frv_nt, std(mfr_nt)];
                sig_idx = strcmp(sig3_fsi(:,1), fsi_labels{iC});
                if any(sig3_fsi{sig_idx,2})
                    sig_fsi_frv_nt = [sig_fsi_frv_nt, std(mfr_nt)];
                 end
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                 % do msn_stuff
                clean_msn = clean_msn + 1;
                tw_mfr = od.msn_res.near_spec{iC}.mfr;
                sig_idx = strcmp(sig3_msn(:,1), msn_labels{iC});
                % Including trials which have min_trial_spikes and have no NaNs
                valid_nothresh = ~isnan(od.msn_res.near_spec{iC}.no_thresh_trialwise_ppc(:,1)); 
                nonan_thresh = ~isnan(od.msn_res.near_spec{iC}.thresh_trialwise_ppc(:,1));
                spike_thresh = od.msn_res.near_spec{iC}.thresh_trialwise_spk_count >= min_trial_spikes;
                valid_thresh = nonan_thresh & spike_thresh;
               
                if sum(valid_thresh) >= min_trials
                    mfr_t = tw_mfr(valid_thresh);
                    mfr_t = (mfr_t - min(mfr_t))/(max(mfr_t) - min(mfr_t)); % Normalize values between 0-1
                    msn_frv = [msn_frv, std(mfr_t)];
                    if any(sig3_msn{sig_idx,2})
                        sig_msn_frv = [sig_msn_frv, std(mfr_t)];
                    end
                end
                
                mfr_nt= tw_mfr(valid_nothresh);
                mfr_nt = (mfr_nt - min(mfr_nt))/(max(mfr_nt) - min(mfr_nt)); % Normalize values between 0-1
                msn_frv_nt = [msn_frv_nt, std(mfr_nt)];
                sig_idx = strcmp(sig3_msn(:,1), msn_labels{iC});
                if any(sig3_msn{sig_idx,2})
                    sig_msn_frv_nt = [sig_msn_frv_nt, std(mfr_nt)];
                 end
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);

%% Plot histograms
fig = figure('WindowState', 'maximized');
ax1 = subplot(2,2,1);
h1 = histogram(fsi_frv_nt, 0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(msn_frv_nt,0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.TickDir = 'out';
ax1.YTick = [0 0.5 1];
ax1.XTick = [0:0.1:0.5];
ax1.XLabel.String = 'Norm Firing Rate SD';
ax1.YLabel.String = 'Proportion';
ax1.XLim = [0 0.5];
ax1.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax1.Box = 'off';
ax1.Title.String = 'No Thesholding, all cells';

ax1 = subplot(2,2,2);
h1 = histogram(fsi_frv, 0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(msn_frv,0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.TickDir = 'out';
ax1.YTick = [0 0.5 1];
ax1.XTick = [0:0.1:0.5];
ax1.XLabel.String = 'Norm Firing Rate SD';
ax1.YLabel.String = 'Proportion';
ax1.XLim = [0 0.5];
ax1.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax1.Box = 'off';
ax1.Title.String = 'Thesholding, all cells';

ax1 = subplot(2,2,3);
h1 = histogram(sig_fsi_frv_nt, 0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(sig_msn_frv_nt,0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.TickDir = 'out';
ax1.YTick = [0 0.5 1];
ax1.XTick = [0:0.1:0.5];
ax1.XLabel.String = 'Norm Firing Rate SD';
ax1.YLabel.String = 'Proportion';
ax1.XLim = [0 0.5];
ax1.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax1.Box = 'off';
ax1.Title.String = 'No Thesholding, Significant cells';

ax1 = subplot(2,2,4);
h1 = histogram(sig_fsi_frv, 0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(sig_msn_frv ,0:0.01:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.TickDir = 'out';
ax1.YTick = [0 0.5 1];
ax1.XTick = [0:0.1:0.5];
ax1.XLabel.String = 'Norm Firing Rate SD';
ax1.YLabel.String = 'Proportion';
ax1.XLim = [0 0.5];
ax1.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax1.Box = 'off';
ax1.Title.String = 'Thesholding, Significant cells';
