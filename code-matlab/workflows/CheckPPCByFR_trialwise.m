%% Calculating and saving the correlation between trialwise firing rates and PPC

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
num_shufs = 1000;

f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
purple = [0.7 0 1]; % Color code for purple
 
fsi_summary = table("label", "delta_corr","delta_sig_pct", ...
    "theta_corr","theta_sig_pct", "beta_corr","beta_sig_pct", ...
    "lo_gamma_corr","lo_gamma_sig_pct", "hi_gamma_corr","hi_gamma_sig_pct");
msn_summary = table("label", "delta_corr","delta_sig_pct", ...
    "theta_corr","theta_sig_pct", "beta_corr","beta_sig_pct", ...
    "lo_gamma_corr","lo_gamma_sig_pct", "hi_gamma_corr","hi_gamma_sig_pct");

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
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ...
                    ~od.fsi_res.near_spec{iC}.flag_no_control_split
                tw_mfr = od.fsi_res.near_spec{iC}.mfr;
                sig_idx = strcmp(sig3_fsi(:,1), fsi_labels{iC});
                this_sig = sig3_fsi{sig_idx,2};
                % Skip this cell if it is not significantly phase locked to any frequency band
                if ~any(this_sig)
                    continue
                end
                % Including trials which have min_trial_spikes and have no NaNs
                nonan_thresh = ~isnan(od.fsi_res.near_spec{iC}.trialwise_unsampled_ppc(:,1));
                spike_thresh = od.fsi_res.near_spec{iC}.thresh_trialwise_spk_count >= min_trial_spikes;
                keep = nonan_thresh & spike_thresh;
                % If the valid number of trials is less than min_trials, then skip
                if sum(keep) < min_trials
                    continue
                end
                clean_fsi = clean_fsi + 1;
                this_row = {fsi_labels{iC}, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
                tw_mfr = tw_mfr(keep);
                tw_ppc = od.fsi_res.near_spec{iC}.trialwise_unsampled_ppc(keep,:);
                this_bands = find(this_sig);
                this_freq = od.fsi_res.near_spec{iC}.freqs;
                this_fppc = zeros(length(this_bands), size(tw_ppc,1));
                tp = length(this_bands);
                fig = figure('WindowState', 'maximized');
                for iB = 1:length(this_bands)
                    f_idx = find(round(this_freq) >= f_list{this_bands(iB)}(1) & ...
                        round(this_freq) <= f_list{this_bands(iB)}(2));
                    this_fppc(iB,:) = mean(tw_ppc(:,f_idx),2)';
                    v_corr = corrcoef(this_fppc(iB,:), tw_mfr);
                    v_corr = v_corr(1,2); % Correlation between PPC and FR
                    % To establish significance of correlation, permute FR w.r.t PPC num_shuf times 
                    % and calculate the correlation for each shuffle
                    this_scorr = zeros(1, num_shufs);
                    for iS = 1:num_shufs
                        s_idx = randperm(length(tw_mfr));
                        scorr = corrcoef(this_fppc(iB,:), tw_mfr(s_idx));
                        this_scorr(iS) = scorr(1,2);
                    end
                    sig_pct = sum(abs(this_scorr) < abs(v_corr))/num_shufs;
                    % Plot trial-wise ppc arranged in increasing order of trial_wise_mfr
                    ax = subplot(tp,1,iB);
                    [~,s_idx] = sort(tw_mfr);
                    plot(tw_mfr(s_idx), this_fppc(iB,s_idx), 'Color', c_list{this_bands(iB)});
                    ax.XAxis.Label.String = 'Trialwise Firing Rate';
                    ax.XAxis.Label.FontSize = 12;
                    ax.YAxis.Label.String = 'PPC';
                    ax.YAxis.Label.FontSize = 12;
                    ax.YAxis.Exponent = 0;
                    ax.Title.String = sprintf("%d Hz - %d Hz", f_list{this_bands(iB)}(1), ...
                        f_list{this_bands(iB)}(2));
                    ax.Title.FontSize = 12;
                    ax.Title.Color = c_list{this_bands(iB)};
                    legend({sprintf('Correlation: %.2f, Sig: %.2f', v_corr, sig_pct)}, ...
                        'Location', 'southwest', 'FontSize', 12);
                    this_row{this_bands(iB)*2} = v_corr;
                    this_row{this_bands(iB)*2+1}  = sig_pct;
                end
                fsi_summary = [fsi_summary; this_row];
                sgtitle(sprintf('FSI: %s', fsi_labels{iC}), 'Interpreter', 'none')
                WriteFig(fig, strcat('D:\RandomVstrAnalysis\temp\FSI_',fsi_labels{iC},'_FRbyPPC'),1)
                close;
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ...
                    ~od.msn_res.near_spec{iC}.flag_no_control_split
                tw_mfr = od.msn_res.near_spec{iC}.mfr;
                sig_idx = strcmp(sig3_msn(:,1), msn_labels{iC});
                this_sig = sig3_msn{sig_idx,2};
                % Skip this cell if it is not significantly phase locked to any frequency band
                if ~any(this_sig)
                    continue
                end
                % Including trials which have min_trial_spikes and have no NaNs
                nonan_thresh = ~isnan(od.msn_res.near_spec{iC}.thresh_trialwise_ppc(:,1));
                spike_thresh = od.msn_res.near_spec{iC}.thresh_trialwise_spk_count >= min_trial_spikes;
                keep = nonan_thresh & spike_thresh;
                % If the valid number of trials is less than min_trials, then skip
                if sum(keep) < min_trials
                    continue
                end
                clean_msn = clean_msn + 1;
                this_row = {msn_labels{iC}, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
                tw_mfr = tw_mfr(keep);
                tw_ppc = od.msn_res.near_spec{iC}.thresh_trialwise_ppc(keep,:);
                this_bands = find(this_sig);
                this_freq = od.msn_res.near_spec{iC}.freqs;
                this_fppc = zeros(length(this_bands), size(tw_ppc,1));
                tp = length(this_bands);
                fig = figure('WindowState', 'maximized');
                for iB = 1:length(this_bands)
                    f_idx = find(round(this_freq) >= f_list{this_bands(iB)}(1) & ...
                        round(this_freq) <= f_list{this_bands(iB)}(2));
                    this_fppc(iB,:) = mean(tw_ppc(:,f_idx),2)';
                    v_corr = corrcoef(this_fppc(iB,:), tw_mfr);
                    v_corr = v_corr(1,2); % Correlation between PPC and FR
                    % To establish significance of correlation, permute FR w.r.t PPC num_shuf times 
                    % and calculate the correlation for each shuffle
                    this_scorr = zeros(1, num_shufs);
                    for iS = 1:num_shufs
                        s_idx = randperm(length(tw_mfr));
                        scorr = corrcoef(this_fppc(iB,:), tw_mfr(s_idx));
                        this_scorr(iS) = scorr(1,2);
                    end
                    sig_pct = sum(abs(this_scorr) < abs(v_corr))/num_shufs;
                    % Plot trial-wise ppc arranged in increasing order of trial_wise_mfr
                    ax = subplot(tp,1,iB);
                    [~,s_idx] = sort(tw_mfr);
                    plot(tw_mfr(s_idx), this_fppc(iB,s_idx), 'Color', c_list{this_bands(iB)});
                    ax.XAxis.Label.String = 'Trialwise Firing Rate';
                    ax.XAxis.Label.FontSize = 12;
                    ax.YAxis.Label.String = 'PPC';
                    ax.YAxis.Label.FontSize = 12;
                    ax.YAxis.Exponent = 0;
                    ax.Title.String = sprintf("%d Hz - %d Hz", f_list{this_bands(iB)}(1), ...
                        f_list{this_bands(iB)}(2));
                    ax.Title.Color = c_list{this_bands(iB)};
                    legend({sprintf('Correlation: %.2f, Sig: %.2f', v_corr, sig_pct)}, ...
                        'Location', 'southwest', 'FontSize', 12);
                    this_row{this_bands(iB)*2} = v_corr;
                    this_row{this_bands(iB)*2+1}  = sig_pct;
                end
                msn_summary = [msn_summary; this_row];
                sgtitle(sprintf('MSN: %s', msn_labels{iC}), 'Interpreter', 'none')
                WriteFig(fig, strcat('D:\RandomVstrAnalysis\temp\MSN_',msn_labels{iC},'_FRbyPPC'),1)
                close;
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
writetable(fsi_summary, strcat('D:\RandomVstrAnalysis\temp\','fsi_summary.xls'));
writetable(msn_summary, strcat('D:\RandomVstrAnalysis\temp\','msn_summary.xls'));
