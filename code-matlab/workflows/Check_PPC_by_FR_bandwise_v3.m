% Script to check if there is an effect of FR by PPC

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
% Load summary data
load('D:\RandomVstrAnalysis\final_results\all_summary.mat');
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
msz = 40; % scatter marker size
f_list = {[3 5], [6 10], [14 25], [40 65], [65 90]};
significant_pairs = cell(1, 5); % Save the score for each of the PPC stuff
p_thresh = 0.99;
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of item
sig3_fsi = {}; % List of significance of each neuron
sig3_msn = {};
hfr_fsi = {};
lfr_fsi = {};
hfr_msn = {};
lfr_msn = {};
ax = {};
fig = figure('WindowState', 'maximized');
for i = 1:2*length(f_list)
    ax{i}  = subplot(2,length(f_list),i);
    hold on;
end    
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

                shuf_ppc = od.fsi_res.near_spec{iC}.shuf_ppc;
                all_ppc = od.fsi_res.near_spec{iC}.subsampled_ppc;
                pct_ppc = sum(shuf_ppc < all_ppc, 1)/size(shuf_ppc,1); % Calculating p_val of PPC for each frequency
                % Test for significance in each of the frequency bands
                this_sig = zeros(1,length(f_list));
                for iF = 1:length(f_list)
                    f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    mean_pct = mean(pct_ppc(f_idx));
                    this_sig(iF) = mean_pct >= p_thresh;
                end

                if any(this_sig) % Go ahead only if any band passes the significance critetion
                    sidx = find(strcmp(fsi_summary.label,fsi_labels{iC}));
                    low_mfr = fsi_summary.lfr_mean(sidx);
                    hi_mfr = fsi_summary.hfr_mean(sidx);
                    lfr_ppc = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc;
                    hfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc;
                    p1_ppc = od.fsi_res.near_p1_spec{iC}.subsampled_ppc;
                    p2_ppc = od.fsi_res.near_p2_spec{iC}.subsampled_ppc;
                    for iF = 1:length(f_list)
                        if ~this_sig(iF)
                            continue
                        end
                        f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                            round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                        plot(ax{iF}, [low_mfr, hi_mfr], [mean(lfr_ppc(f_idx)), mean(hfr_ppc(f_idx))], 'color', c_list{iF});
                        
                        pct_hfr_ppc = sum(shuf_ppc < hfr_ppc, 1)/size(shuf_ppc,1);
                        mean_hfr_pct = mean(pct_hfr_ppc(f_idx));
%                         if mean_hfr_pct >= p_thresh
%                            scatter(ax{iF}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF})
%                         end
% 
%                         pct_lfr_ppc = sum(shuf_ppc < lfr_ppc, 1)/size(shuf_ppc,1);
%                         mean_lfr_pct = mean(pct_lfr_ppc(f_idx));
%                         if mean_lfr_pct >= p_thresh
%                            scatter(ax{iF}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF})
%                         end

                        % If the difference is significant, make dark circles
                        this_diff = abs(mean(hfr_ppc(f_idx)) - mean(lfr_ppc(f_idx)));
                        this_cdiff = abs(mean(p1_ppc(:,f_idx),2) - mean(p2_ppc(:,f_idx),2));
                        if this_diff >= prctile(this_cdiff,100*p_thresh)
                            scatter(ax{iF}, 0.5*(low_mfr+hi_mfr), ...
                                0.5*(mean(lfr_ppc(f_idx))+mean(hfr_ppc(f_idx))), ...
                                msz, c_list{iF}, 'Filled')
                        end
                    end 
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

                shuf_ppc = od.msn_res.near_spec{iC}.shuf_ppc;
                all_ppc = od.msn_res.near_spec{iC}.ppc';
                pct_ppc = sum(shuf_ppc < all_ppc, 1)/size(shuf_ppc,1); % Calculating p_val of PPC for each frequency
                % Test for significance in each of the frequency bands
                this_sig = zeros(1,length(f_list));
                for iF = 1:length(f_list)
                    f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    mean_pct = mean(pct_ppc(f_idx));
                    this_sig(iF) = mean_pct >= p_thresh;
                end

                if any(this_sig) % Go ahead only if any band passes the significance critetion
                    sidx = find(strcmp(msn_summary.label,msn_labels{iC}));
                    low_mfr = msn_summary.lfr_mean(sidx);
                    hi_mfr = msn_summary.hfr_mean(sidx);
                    lfr_ppc = od.msn_res.near_lfr_spec{iC}.ppc';
                    hfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc';
                    p1_ppc = od.msn_res.near_p1_spec{iC}.ppc;
                    p2_ppc = od.msn_res.near_p2_spec{iC}.ppc;
                    for iF = 1:length(f_list)
                        if ~this_sig(iF)
                            continue
                        end
                        f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                            round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                        plot(ax{iF}, [low_mfr, hi_mfr], [mean(lfr_ppc(f_idx)), mean(hfr_ppc(f_idx))], 'color', c_list{iF});

                        pct_hfr_ppc = sum(shuf_ppc < hfr_ppc, 1)/size(shuf_ppc,1); % Calculating p_val of PPC for each frequency
                        mean_hfr_pct = mean(pct_hfr_ppc(f_idx));
%                         if mean_hfr_pct >= p_thresh
%                            scatter(ax{iF}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF})
%                         end
% 
%                         pct_lfr_ppc = sum(shuf_ppc < lfr_ppc, 1)/size(shuf_ppc,1); % Calculating p_val of PPC for each frequency
%                         mean_lfr_pct = mean(pct_lfr_ppc(f_idx));
%                         if mean_lfr_pct >= p_thresh
%                            scatter(ax{iF}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF})
%                         end

                        % If the difference is significant, make dark circles
                        this_diff = abs(mean(hfr_ppc(f_idx)) - mean(lfr_ppc(f_idx)));
                        this_cdiff = abs(mean(p1_ppc(:,f_idx),2) - mean(p2_ppc(:,f_idx),2));
                        if this_diff >= prctile(this_cdiff,100*p_thresh)
                            scatter(ax{iF}, 0.5*(low_mfr+hi_mfr), ...
                                0.5*(mean(lfr_ppc(f_idx))+mean(hfr_ppc(f_idx))), ...
                                msz, c_list{iF}, 'Filled')
                        end
                    end 
                end
            end
        end
    end
end
% Set Axis Labels and names
for iF=1:length(f_list)
    ax{iF}.Title.String = sprintf("FSI %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2));
    ax{iF}.Title.Color =  c_list{iF};
    ax{iF}.YAxis.Exponent = 0;
    ax{iF}.XLim = [0 80];
    ax{iF+length(f_list)}.Title.String = sprintf("MSN %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2));
    ax{iF+length(f_list)}.Title.Color =  c_list{iF};
    ax{iF+length(f_list)}.YAxis.Exponent = 0;
    ax{iF+length(f_list)}.XLim = [0 12];
end

fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
%% find multi-band phase-locked MSNs and FSIs
% load('D:\RandomVstrAnalysis\final_results\significantly_phase_locked.mat');
% qm = cell2mat(sig3_msn(:,2));
% qf = cell2mat(sig3_fsi(:,2));
% sm = sum(qm,2);
% sf = sum(qf,2);
% plf = find(sf > 1);
% plm = find(sm > 1);
%% Check how many of these neurons also show significant firing rate based modulation
% load('D:\RandomVstrAnalysis\final_results\significance_v2.mat');
% % percentile based significance
% pct = cell2mat(sig2_fsi(:,2));
% sig_pct = pct(plf,:);
% sig_band = zeros(length(plf), length(f_list));
% for iS = 1:length(sig_pct)
%     for iF = 1:length(f_list)
%         f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
%             round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
%         if 
%     end
% end
%% 

