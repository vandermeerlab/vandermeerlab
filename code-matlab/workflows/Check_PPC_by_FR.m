% Script to check if there is an effect of FR by PPC

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
% load significant cells
load('significance.mat'); % Label, HFR PPC high, LFR PPC high, HFR STS high, LFR STS high
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
msz = 40; % scatter marker size
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
significant_pairs = cell(1, 5); % Save the score for each of the PPC stuff
p_thresh = 0.95;
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of item
ax = {};
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
                % Load the PPC_by_FR files to extract the shuffle_ppc
                shuf_data = load(strcat('D:\RandomVstrAnalysis\PPCbyFR\Old\', fsi_labels{iC}, '_FSI_ppcbyfr.mat'));
                shuf_ppc = shuf_data.od.shuf_ppc;
                all_ppc = od.fsi_res.near_spec{iC}.subsampled_ppc;
                pctile_ppc = sum(shuf_ppc < all_ppc, 1)/shuf_data.od.nshufs; % Calculating p_val of PPC for each frequency
                % Test for significance in each of the frequency bands
                this_sig = zeros(1,length(f_list));
                for iF = 1:length(f_list)
                    f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    mean_pctile = mean(pctile_ppc(f_idx));
                    this_sig(iF) = mean_pctile >= p_thresh;
                end
                if any(this_sig) % Go ahead only if any significant value is present
                    % Calculate the mfr
                    mfr = od.fsi_res.near_spec{iC}.mfr;
                    nz_mfr = mfr(mfr > 0);
                    low_mfr = mean(nz_mfr((nz_mfr <= od.fsi_res.near_spec{iC}.fr_thresh)));
                    hi_mfr = mean(nz_mfr(nz_mfr > od.fsi_res.near_spec{iC}.fr_thresh));
                    lfr_ppc = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc;
                    hfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc;
                    for iF = 1:length(f_list)
                        if ~this_sig(iF)
                            continue
                        end
                        f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                        plot(ax{iF}, [low_mfr, hi_mfr], [mean(lfr_ppc(f_idx)), mean(hfr_ppc(f_idx))], 'color', c_list{iF});
                        % Check if this was a significant cell
                        sig_idx = strcmp(sig_fsi(:,1),fsi_labels{iC});
                        if sig_fsi{sig_idx, 2} ~= 0
                            if (round(sig_fsi{sig_idx, 2}) >= f_list{iF}(1)) && (round(sig_fsi{sig_idx, 2}) <= f_list{iF}(2))
                                scatter(ax{iF}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF}, 'filled')
                            else
                                scatter(ax{iF}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF})
                            end
                        end
                        if sig_fsi{sig_idx, 3} ~= 0
                            if (round(sig_fsi{sig_idx, 3}) >= f_list{iF}(1)) && (round(sig_fsi{sig_idx, 3}) <= f_list{iF}(2))
                                scatter(ax{iF}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF}, 'filled')
                            else
                                scatter(ax{iF}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF})
                            end
                        end
                        if sig_fsi{sig_idx, 4} ~= 0
                            if (round(sig_fsi{sig_idx, 2}) >= f_list{iF}(1)) && (round(sig_fsi{sig_idx, 2}) <= f_list{iF}(2))
                                scatter(ax{iF}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF}, 'filled', 'd')
                            else
                                scatter(ax{iF}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF}, 'd')
                            end
                        end
                        if sig_fsi{sig_idx, 5} ~= 0
                            if (round(sig_fsi{sig_idx, 3}) >= f_list{iF}(1)) && (round(sig_fsi{sig_idx, 3}) <= f_list{iF}(2))
                                scatter(ax{iF}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF}, 'filled', 'd')
                            else
                                scatter(ax{iF}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF}, 'd')
                            end
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
                clean_msn = clean_fsi + 1;
                % Load the PPC_by_FR files to extract the shuffle_ppc
                shuf_data = load(strcat('D:\RandomVstrAnalysis\PPCbyFR\Old\', msn_labels{iC}, '_MSN_ppcbyfr.mat'));
                shuf_ppc = shuf_data.od.shuf_ppc;
                all_ppc = od.msn_res.near_spec{iC}.ppc';
                pctile_ppc = sum(shuf_ppc < all_ppc, 1)/shuf_data.od.nshufs; % Calculating p_val of PPC for each frequency
                % Test for significance in each of the frequency bands
                this_sig = zeros(1,length(f_list));
                for iF = 1:length(f_list)
                    f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    mean_pctile = mean(pctile_ppc(f_idx));
                    this_sig(iF) = mean_pctile >= p_thresh;
                end
                if any(this_sig) % Go ahead only if any significant value is present
                    % Calculate the mfr
                    mfr = od.msn_res.near_spec{iC}.mfr;
                    nz_mfr = mfr(mfr > 0);
                    low_mfr = mean(nz_mfr((nz_mfr <= od.msn_res.near_spec{iC}.fr_thresh)));
                    hi_mfr = mean(nz_mfr(nz_mfr > od.msn_res.near_spec{iC}.fr_thresh));
                    lfr_ppc = od.msn_res.near_lfr_spec{iC}.ppc';
                    hfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc';
                    for iF = 1:length(f_list)
                        if ~this_sig(iF)
                            continue
                        end
                        f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                        plot(ax{iF+length(f_list)}, [low_mfr, hi_mfr], [mean(lfr_ppc(f_idx)), mean(hfr_ppc(f_idx))], 'color', c_list{iF});
                        % Check if this was a significant cell
                        sig_idx = strcmp(sig_msn(:,1),msn_labels{iC});
                        if sig_msn{sig_idx, 2} ~= 0
                            if (round(sig_msn{sig_idx, 2}) >= f_list{iF}(1)) && (round(sig_msn{sig_idx, 2}) <= f_list{iF}(2))
                                scatter(ax{iF+length(f_list)}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF}, 'filled')
                            else
                                scatter(ax{iF+length(f_list)}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF})
                            end
                        end
                        if sig_msn{sig_idx, 3} ~= 0
                            if (round(sig_msn{sig_idx, 3}) >= f_list{iF}(1)) && (round(sig_msn{sig_idx, 3}) <= f_list{iF}(2))
                                scatter(ax{iF+length(f_list)}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF}, 'filled')
                            else
                                scatter(ax{iF+length(f_list)}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF})
                            end
                        end
                        if sig_msn{sig_idx, 4} ~= 0
                            if (round(sig_msn{sig_idx, 2}) >= f_list{iF}(1)) && (round(sig_msn{sig_idx, 2}) <= f_list{iF}(2))
                                scatter(ax{iF+length(f_list)}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF}, 'filled', 'd')
                            else
                                scatter(ax{iF+length(f_list)}, hi_mfr, mean(hfr_ppc(f_idx)), msz, c_list{iF}, 'd')
                            end
                        end
                        if sig_msn{sig_idx, 5} ~= 0
                            if (round(sig_msn{sig_idx, 3}) >= f_list{iF}(1)) && (round(sig_msn{sig_idx, 3}) <= f_list{iF}(2))
                                scatter(ax{iF+length(f_list)}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF}, 'filled', 'd')
                            else
                                scatter(ax{iF+length(f_list)}, low_mfr, mean(lfr_ppc(f_idx)), msz, c_list{iF}, 'd')
                            end
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

