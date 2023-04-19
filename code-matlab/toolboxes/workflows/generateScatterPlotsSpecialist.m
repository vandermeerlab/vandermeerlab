% Generates Scatter Plots of Difference in Peak Frequencies vs difference
% in mean firing rates

%SD is not being calculated properly
fsi_near_lfr_ppc_peaks = [];
fsi_near_hfr_ppc_peaks = [];
fsi_near_hfr_pos_peak_ptile = [];
fsi_near_hfr_pos_peak_ratio = [];
fsi_near_lfr_pos_peak_ptile = [];
fsi_near_lfr_pos_peak_ratio = [];
fsi_near_hfr_neg_peak_ptile = [];
fsi_near_hfr_neg_peak_ratio = [];
fsi_near_lfr_neg_peak_ptile = [];
fsi_near_lfr_neg_peak_ratio = [];
fsi_near_lfr_frs = [];
fsi_near_hfr_frs = [];
fsi_max_difs = [];
fsi_min_difs = [];
fsi_near_ppc_peaks = [];
fsi_near_ppc_maxval = [];
fsi_near_hfr_ppc_maxval = [];
fsi_near_lfr_ppc_maxval = [];
valid_fsi_labels = [];
sig_dif_fsi = []; %0 if none, 1 if lfr>hfr, 2 if hfr>lfr, 3 if both

msn_near_lfr_ppc_peaks = [];
msn_near_hfr_ppc_peaks = [];
msn_near_hfr_pos_peak_ptile = [];
msn_near_hfr_pos_peak_ratio = [];
msn_near_lfr_pos_peak_ptile = [];
msn_near_lfr_pos_peak_ratio = [];
msn_near_hfr_neg_peak_ptile = [];
msn_near_hfr_neg_peak_ratio = [];
msn_near_lfr_neg_peak_ptile = [];
msn_near_lfr_neg_peak_ratio = [];
msn_near_lfr_frs = [];
msn_near_hfr_frs = [];
msn_max_difs = [];
msn_min_difs = [];
msn_near_ppc_peaks = [];
msn_near_ppc_maxval = [];
msn_near_hfr_ppc_maxval = [];
msn_near_lfr_ppc_maxval = [];
valid_msn_labels = [];
sig_dif_msn = []; %0 if none, 1 if lfr>hfr, 2 if hfr>lfr, 3 if both

cd('D:\RandomVstrAnalysis\ft_results');

rats = {'R117','R119','R131','R132'};
fw = [5,100];
pk_thresh = -1;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        fsi_labels  = od.label(od.cell_type == 2);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC},'flag_no_control_split') && ...
                        ~od.fsi_res.near_spec{iC}.flag_no_control_split  
                
                flag_near_lfr_ppc_peak = true;
                flag_near_hfr_ppc_peak = true;
                flag_near_ppc_peak = true;
                
                %Peaks in near_spec
                lf = find(od.fsi_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_spec{iC}.subsampled_ppc(lf:rf));
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_ppc_val = max_val;
                        near_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_ppc_peak = false;
                    end
                else
                    flag_near_ppc_peak = false; 
                end   

                % Peaks in near_lfr_spec
                lf = find(od.fsi_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf:rf));
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_ppc_val = max_val;
                        near_lfr_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_ppc_peak = false;
                    end
                else
                    flag_near_lfr_ppc_peak = false; 
                end
                
                % Peaks in near_hfr_spec               
                lf = find(od.fsi_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf:rf));
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_ppc_val = max_val;
                        near_hfr_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_ppc_peak = false;
                    end
                else
                    flag_near_hfr_ppc_peak = false; 
                end
           
                all_ppc_peaks = flag_near_hfr_ppc_peak && ...
                    flag_near_lfr_ppc_peak;

                % Skip plotting only if no peaks  
                if ~all_ppc_peaks  
                    close all;
                    continue;
                end
                
                % Calculate mean mfr and range
                nz_trials = od.fsi_res.near_spec{iC}.mfr > 0;
                lfr_trials = od.fsi_res.near_spec{iC}.mfr <= od.fsi_res.near_spec{iC}.fr_thresh;
                lfr_trials = lfr_trials & nz_trials;
                hfr_trials = ~lfr_trials;
                hfr_trials = hfr_trials & nz_trials;
                near_lfr_spec_mean = sum(od.fsi_res.near_spec{iC}.trial_spk_count(lfr_trials))/sum(od.fsi_res.near_spec{iC}.trial_spk_count(lfr_trials)'./od.fsi_res.near_spec{iC}.mfr(lfr_trials));
                near_lfr_spec_range = [min(od.fsi_res.near_spec{iC}.mfr(lfr_trials)), max(od.fsi_res.near_spec{iC}.mfr(lfr_trials))];
                near_hfr_spec_mean = sum(od.fsi_res.near_spec{iC}.trial_spk_count(hfr_trials))/sum(od.fsi_res.near_spec{iC}.trial_spk_count(hfr_trials)'./od.fsi_res.near_spec{iC}.mfr(hfr_trials));
                near_hfr_spec_range = [min(od.fsi_res.near_spec{iC}.mfr(hfr_trials)), max(od.fsi_res.near_spec{iC}.mfr(hfr_trials))];
                near_spec_mean = sum(od.fsi_res.near_spec{iC}.trial_spk_count(nz_trials))/sum(od.fsi_res.near_spec{iC}.trial_spk_count(nz_trials)'./od.fsi_res.near_spec{iC}.mfr(nz_trials));
                near_spec_range = [min(od.fsi_res.near_spec{iC}.mfr(nz_trials)), max(od.fsi_res.near_spec{iC}.mfr(nz_trials))];
                
                %If both HFR and LFR dif have a positive and negative diff
                lf = find(od.fsi_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                hfr_lfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf:rf) - ...
                    od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf:rf);
                p1_p2_ppc = od.fsi_res.near_p1_spec{iC}.subsampled_ppc(:,lf:rf) - ...
                    od.fsi_res.near_p2_spec{iC}.subsampled_ppc(:,lf:rf);
                mean_ppc = mean(p1_p2_ppc, 1);
                sd_ppc = std(p1_p2_ppc, 1);
                pos_difs = hfr_lfr_ppc - (mean_ppc + 3*sd_ppc);
                neg_difs =(mean_ppc - 3*sd_ppc) - hfr_lfr_ppc;
                [dif_max, max_loc] = max(pos_difs);
                [dif_min, min_loc] = max(neg_difs);
                
                %Save the percentile and the ratio of the highest and
                %lowest ppc_dif
                this_hfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf:rf);
                this_lfr_ppc = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf:rf);
                pos_sig_hfr_ptile = sum(this_hfr_ppc(max_loc) >= this_hfr_ppc)/length(this_hfr_ppc);
                pos_sig_lfr_ptile = sum(this_lfr_ppc(max_loc) >= this_lfr_ppc)/length(this_lfr_ppc);
                pos_sig_hfr_ratio = this_hfr_ppc(max_loc)/max(this_hfr_ppc);
                pos_sig_lfr_ratio = this_lfr_ppc(max_loc)/max(this_lfr_ppc);
                
                neg_sig_hfr_ptile = sum(this_hfr_ppc(min_loc) >= this_hfr_ppc)/length(this_hfr_ppc);
                neg_sig_lfr_ptile = sum(this_lfr_ppc(min_loc) >= this_lfr_ppc)/length(this_lfr_ppc);
                neg_sig_hfr_ratio = this_hfr_ppc(min_loc)/max(this_hfr_ppc);
                neg_sig_lfr_ratio = this_lfr_ppc(min_loc)/max(this_lfr_ppc);
         
                fsi_near_hfr_pos_peak_ptile = [fsi_near_hfr_pos_peak_ptile pos_sig_hfr_ptile];
                fsi_near_hfr_pos_peak_ratio = [fsi_near_hfr_pos_peak_ratio pos_sig_hfr_ratio];
                fsi_near_lfr_pos_peak_ptile = [fsi_near_lfr_pos_peak_ptile pos_sig_lfr_ptile];
                fsi_near_lfr_pos_peak_ratio = [fsi_near_lfr_pos_peak_ratio pos_sig_lfr_ratio];
                
                fsi_near_hfr_neg_peak_ptile = [fsi_near_hfr_neg_peak_ptile neg_sig_hfr_ptile];
                fsi_near_hfr_neg_peak_ratio = [fsi_near_hfr_neg_peak_ratio neg_sig_hfr_ratio];
                fsi_near_lfr_neg_peak_ptile = [fsi_near_lfr_neg_peak_ptile neg_sig_lfr_ptile];
                fsi_near_lfr_neg_peak_ratio = [fsi_near_lfr_neg_peak_ratio neg_sig_lfr_ratio];
                
                %Mark significance of ppc_dif
                if sum(hfr_lfr_ppc >= mean_ppc + 3*sd_ppc) ~= 0 & ...
                        sum(hfr_lfr_ppc <= mean_ppc - 3*sd_ppc) ~= 0
                    sig_dif_fsi = [sig_dif_fsi 3];
                elseif sum(hfr_lfr_ppc >= mean_ppc + 3*sd_ppc) ~= 0
                    sig_dif_fsi = [sig_dif_fsi 2];
                elseif sum(hfr_lfr_ppc <= mean_ppc - 3*sd_ppc) ~= 0
                    sig_dif_fsi = [sig_dif_fsi 1];
                else
                    sig_dif_fsi = [sig_dif_fsi 0];
                end
                
                % if peaks exist 
                if all_ppc_peaks
                    near_lfr_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(fw(1)+near_lfr_ppc_pk-1);
                    near_hfr_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(fw(1)+near_hfr_ppc_pk-1);
                    
                    fsi_near_ppc_peaks = [fsi_near_ppc_peaks near_ppc_pk];    
                    fsi_near_lfr_ppc_peaks = [fsi_near_lfr_ppc_peaks, near_lfr_ppc_pk];
                    fsi_near_hfr_ppc_peaks = [fsi_near_hfr_ppc_peaks, near_hfr_ppc_pk];
                    
                    fsi_near_ppc_maxval = [fsi_near_ppc_maxval near_ppc_val];
                    fsi_near_lfr_ppc_maxval = [fsi_near_lfr_ppc_maxval near_lfr_ppc_val];
                    fsi_near_hfr_ppc_maxval = [fsi_near_hfr_ppc_maxval near_hfr_ppc_val];
                    
                    fsi_near_lfr_frs = [fsi_near_lfr_frs, near_lfr_spec_mean];
                    fsi_near_hfr_frs = [fsi_near_hfr_frs, near_hfr_spec_mean];
                    valid_fsi_labels = [valid_fsi_labels fsi_labels(iC)];
                    
                    fsi_max_difs = [fsi_max_difs od.fsi_res.near_hfr_spec{iC}.freqs(fw(1)+max_loc-1)];
                    fsi_min_difs = [fsi_min_difs od.fsi_res.near_lfr_spec{iC}.freqs(fw(1)+min_loc-1)];
                end
            else
                continue
            end 
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC},'flag_no_control_split') && ...
                    ~od.msn_res.near_spec{iC}.flag_no_control_split  
                
                flag_near_lfr_ppc_peak = true;
                flag_near_hfr_ppc_peak = true;
                
                %Peaks in near_spec
                lf = find(od.msn_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_spec{iC}.ppc(lf:rf));
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_ppc_val = max_val;
                        near_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_ppc_peak = false;
                    end
                else
                    flag_near_ppc_peak = false; 
                end   
                
                % Peaks in near_lfr_spec              
                lf = find(od.msn_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_lfr_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_lfr_spec{iC}.ppc(lf:rf));
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_ppc_val = max_val;
                        near_lfr_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_ppc_peak = false;
                    end
                else
                    flag_near_lfr_ppc_peak = false; 
                end
                
                % Peaks in near_hfr_spec              
                lf = find(od.msn_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_hfr_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_hfr_spec{iC}.ppc(lf:rf));
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_ppc_val = max_val;
                        near_hfr_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_ppc_peak = false;
                    end
                else
                    flag_near_hfr_ppc_peak = false; 
                end
                
                all_ppc_peaks = flag_near_hfr_ppc_peak && ...
                    flag_near_lfr_ppc_peak;

                % Skip plotting only if no peaks  
                if ~all_ppc_peaks  
                    close all;
                    continue;
                end
                
                % Calculate mean mfr and range
                nz_trials = od.msn_res.near_spec{iC}.mfr > 0;
                lfr_trials = od.msn_res.near_spec{iC}.mfr <= od.msn_res.near_spec{iC}.fr_thresh;
                lfr_trials = lfr_trials & nz_trials;
                hfr_trials = ~lfr_trials;
                hfr_trials = hfr_trials & nz_trials;
                near_lfr_spec_mean = sum(od.msn_res.near_spec{iC}.trial_spk_count(lfr_trials))/sum(od.msn_res.near_spec{iC}.trial_spk_count(lfr_trials)'./od.msn_res.near_spec{iC}.mfr(lfr_trials));
                near_lfr_spec_range = [min(od.msn_res.near_spec{iC}.mfr(lfr_trials)), max(od.msn_res.near_spec{iC}.mfr(lfr_trials))];
                near_hfr_spec_mean = sum(od.msn_res.near_spec{iC}.trial_spk_count(hfr_trials))/sum(od.msn_res.near_spec{iC}.trial_spk_count(hfr_trials)'./od.msn_res.near_spec{iC}.mfr(hfr_trials));
                near_hfr_spec_range = [min(od.msn_res.near_spec{iC}.mfr(hfr_trials)), max(od.msn_res.near_spec{iC}.mfr(hfr_trials))];
                
                %If both HFR and LFR dif have a positive and negative diff
                lf = find(od.msn_res.near_spec{iC}.freqs >= fw(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= fw(2), 1, 'last');
                hfr_lfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc(lf:rf) - ...
                    od.msn_res.near_lfr_spec{iC}.ppc(lf:rf);
                %Need to transpose this to convert to suitable dimensions
                hfr_lfr_ppc = hfr_lfr_ppc';
                p1_p2_ppc = od.msn_res.near_p1_spec{iC}.ppc(:,lf:rf) - ...
                    od.msn_res.near_p2_spec{iC}.ppc(:,lf:rf);
                mean_ppc = mean(p1_p2_ppc, 1);
                sd_ppc = std(p1_p2_ppc, 1);
                pos_difs = hfr_lfr_ppc - (mean_ppc + 3*sd_ppc);
                neg_difs =(mean_ppc - 3*sd_ppc) - hfr_lfr_ppc;
                [dif_max, max_loc] = max(pos_difs);
                [dif_min, min_loc] = max(neg_difs);
                
                %Save the percentile and the ratio of the highest and
                %lowest ppc_dif
                this_hfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc(lf:rf);
                this_lfr_ppc = od.msn_res.near_lfr_spec{iC}.ppc(lf:rf);
                pos_sig_hfr_ptile = sum(this_hfr_ppc(max_loc) >= this_hfr_ppc)/length(this_hfr_ppc);
                pos_sig_lfr_ptile = sum(this_lfr_ppc(max_loc) >= this_lfr_ppc)/length(this_lfr_ppc);
                pos_sig_hfr_ratio = this_hfr_ppc(max_loc)/max(this_hfr_ppc);
                pos_sig_lfr_ratio = this_lfr_ppc(max_loc)/max(this_lfr_ppc);
                
                neg_sig_hfr_ptile = sum(this_hfr_ppc(min_loc) >= this_hfr_ppc)/length(this_hfr_ppc);
                neg_sig_lfr_ptile = sum(this_lfr_ppc(min_loc) >= this_lfr_ppc)/length(this_lfr_ppc);
                neg_sig_hfr_ratio = this_hfr_ppc(min_loc)/max(this_hfr_ppc);
                neg_sig_lfr_ratio = this_lfr_ppc(min_loc)/max(this_lfr_ppc);
         
                msn_near_hfr_pos_peak_ptile = [msn_near_hfr_pos_peak_ptile pos_sig_hfr_ptile];
                msn_near_hfr_pos_peak_ratio = [msn_near_hfr_pos_peak_ratio pos_sig_hfr_ratio];
                msn_near_lfr_pos_peak_ptile = [msn_near_lfr_pos_peak_ptile pos_sig_lfr_ptile];
                msn_near_lfr_pos_peak_ratio = [msn_near_lfr_pos_peak_ratio pos_sig_lfr_ratio];
                
                msn_near_hfr_neg_peak_ptile = [msn_near_hfr_neg_peak_ptile neg_sig_hfr_ptile];
                msn_near_hfr_neg_peak_ratio = [msn_near_hfr_neg_peak_ratio neg_sig_hfr_ratio];
                msn_near_lfr_neg_peak_ptile = [msn_near_lfr_neg_peak_ptile neg_sig_lfr_ptile];
                msn_near_lfr_neg_peak_ratio = [msn_near_lfr_neg_peak_ratio neg_sig_lfr_ratio];
                            
                %Mark significance of ppc_dif
                if sum(hfr_lfr_ppc >= mean_ppc + 3*sd_ppc) ~= 0 & ...
                        sum(hfr_lfr_ppc <= mean_ppc - 3*sd_ppc) ~= 0
                    sig_dif_msn = [sig_dif_msn 3];
                elseif sum(hfr_lfr_ppc >= mean_ppc + 3*sd_ppc) ~= 0
                    sig_dif_msn = [sig_dif_msn 2];
                elseif sum(hfr_lfr_ppc <= mean_ppc - 3*sd_ppc) ~= 0
                    sig_dif_msn = [sig_dif_msn 1];
                else
                    sig_dif_msn = [sig_dif_msn 0];
                end
                
                % if peaks exist
                if all_ppc_peaks
                    near_lfr_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(fw(1)+near_lfr_ppc_pk-1);
                    near_hfr_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(fw(1)+near_hfr_ppc_pk-1);  

                    msn_near_ppc_peaks = [msn_near_ppc_peaks near_ppc_pk];    
                    msn_near_lfr_ppc_peaks = [msn_near_lfr_ppc_peaks, near_lfr_ppc_pk];
                    msn_near_hfr_ppc_peaks = [msn_near_hfr_ppc_peaks, near_hfr_ppc_pk];
                    
                    msn_near_ppc_maxval = [msn_near_ppc_maxval near_ppc_val];
                    msn_near_lfr_ppc_maxval = [msn_near_lfr_ppc_maxval near_lfr_ppc_val];
                    msn_near_hfr_ppc_maxval = [msn_near_hfr_ppc_maxval near_hfr_ppc_val];
                    
                    msn_near_lfr_frs = [msn_near_lfr_frs, near_lfr_spec_mean];
                    msn_near_hfr_frs = [msn_near_hfr_frs, near_hfr_spec_mean];
                    valid_msn_labels = [valid_msn_labels msn_labels(iC)];
                    
                    msn_max_difs = [msn_max_difs od.msn_res.near_hfr_spec{iC}.freqs(fw(1)+max_loc-1)];
                    msn_min_difs = [msn_min_difs od.msn_res.near_lfr_spec{iC}.freqs(fw(1)+min_loc-1)];
                end
            else
                continue
            end
        end
    end
end

%% Plot Distance betweem Peak Frequencies vs Difference in Mean Firing rates in MSNs and FSIs
fig_1 = figure('WindowState', 'maximized');
s1 = scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_ppc_peaks - msn_near_lfr_ppc_peaks);
s1.Marker = 'o';
s1.MarkerFaceColor = [0.8500 0.3250 0.0980];
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeColor = [0.8500 0.3250 0.0980];
s1.MarkerEdgeAlpha = 0.5;
s1.SizeData = 150;
hold on;
s2 = scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_near_hfr_ppc_peaks - fsi_near_lfr_ppc_peaks);
s2.Marker = 'o';
s2.MarkerFaceColor = 'blue';
s2.MarkerFaceAlpha = 0.5;
s2.MarkerEdgeColor = 'blue';
s2.MarkerEdgeAlpha = 0.5;
s2.SizeData = 150;
ax1 = gca(fig_1);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.XLim = [0,10];
ax1.YLim = [-120 120];
ax1.YLabel.String = '\Delta Peak Frequency';
ax1.XLabel.String = '\Delta Firing Rate';
ax1.XLabel.FontSize = 24;
ax1.YLabel.FontSize = 24;
leg = legend({'MSN', 'FSI'});
leg.FontName = 'Arial';
leg.FontSize = 20;
leg.FontWeight = 'bold';
box off;
%% Plot scatter with FSI labels
fig_2 = figure('WindowState', 'maximized');
hold on;
for iF = 1:length(valid_fsi_labels)
    text(fsi_near_hfr_frs(iF) - fsi_near_lfr_frs(iF), ...
        fsi_near_hfr_ppc_peaks(iF) - fsi_near_lfr_ppc_peaks(iF), ...
        valid_fsi_labels(iF), 'Interpreter', 'none');
end
xlim([0,10])
ylim([-100 100])
title('FSI Label scatter');
%% Plot scatter with MSN labels
fig_3 = figure('WindowState', 'maximized');
hold on;
for iF = 1:length(valid_msn_labels)
    text(msn_near_hfr_frs(iF) - msn_near_lfr_frs(iF), ...
        msn_near_hfr_ppc_peaks(iF) - msn_near_lfr_ppc_peaks(iF), ...
        valid_msn_labels(iF), 'Interpreter', 'none');
end
xlim([0,10])
ylim([-100 100])
title('MSN Label scatter');
%% Plot Distance betweem Peak Frequencies vs Difference in Mean Firing rates in MSNs and FSIs
fig_4 = figure('WindowState', 'maximized');
s1 = scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_max_difs - msn_min_difs);
s1.Marker = 'o';
s1.MarkerFaceColor = [0.8500 0.3250 0.0980];
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeColor = [0.8500 0.3250 0.0980];
s1.MarkerEdgeAlpha = 0.5;
s1.SizeData = 150;
hold on;
s2 = scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_max_difs - fsi_min_difs);
s2.Marker = 'o';
s2.MarkerFaceColor = 'blue';
s2.MarkerFaceAlpha = 0.5;
s2.MarkerEdgeColor = 'blue';
s2.MarkerEdgeAlpha = 0.5;
s2.SizeData = 150;
ax1 = gca(fig_1);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.XLim = [0,10];
ax1.YLim = [-120 120];
ax1.YLabel.String = '\Delta Peak Frequency';
ax1.XLabel.String = '\Delta Firing Rate';
ax1.XLabel.FontSize = 24;
ax1.YLabel.FontSize = 24;
leg = legend({'MSN', 'FSI'});
leg.FontName = 'Arial';
leg.FontSize = 20;
leg.FontWeight = 'bold';
box off;

%% Plot scatters of different methods of peak differences
fig_4 = figure('WindowState', 'maximized');
s1 = scatter(msn_near_hfr_ppc_peaks - msn_near_lfr_ppc_peaks, ...
    msn_max_difs - msn_min_difs);
keyboard;
s2 = scatter(fsi_near_hfr_ppc_peaks - fsi_near_lfr_ppc_peaks, ...
    fsi_max_difs - fsi_min_difs);

%% Significant difference sandbox
% 0 if none, 1 if lfr>hfr, 2 if hfr>lfr, 3 if both
sig_fsi_labels = valid_fsi_labels(sig_dif_fsi ~= 0);
sig_msn_labels = valid_msn_labels(sig_dif_msn ~= 0);
sig_fsi_type = sig_dif_fsi(sig_dif_fsi~= 0);
sig_msn_type = sig_dif_msn(sig_dif_msn~= 0);
sig_msn_max_difs = msn_max_difs(sig_dif_msn~= 0);
sig_msn_min_difs = msn_min_difs(sig_dif_msn~= 0);
sig_fsi_max_difs = fsi_max_difs(sig_dif_fsi~= 0);
sig_fsi_min_difs = fsi_min_difs(sig_dif_fsi~= 0);
sig_fsi_pos_hfr_peak_ptile = fsi_near_hfr_pos_peak_ptile(sig_dif_fsi~= 0);
sig_fsi_pos_lfr_peak_ptile = fsi_near_lfr_pos_peak_ptile(sig_dif_fsi~= 0);
sig_fsi_pos_hfr_peak_ratio = fsi_near_hfr_pos_peak_ratio(sig_dif_fsi~= 0);
sig_fsi_pos_lfr_peak_ratio = fsi_near_lfr_pos_peak_ratio(sig_dif_fsi~= 0);
sig_fsi_neg_hfr_peak_ptile = fsi_near_hfr_neg_peak_ptile(sig_dif_fsi~= 0);
sig_fsi_neg_lfr_peak_ptile = fsi_near_lfr_neg_peak_ptile(sig_dif_fsi~= 0);
sig_fsi_neg_hfr_peak_ratio = fsi_near_hfr_neg_peak_ratio(sig_dif_fsi~= 0);
sig_fsi_neg_lfr_peak_ratio = fsi_near_lfr_neg_peak_ratio(sig_dif_fsi~= 0);
sig_msn_pos_hfr_peak_ptile = msn_near_hfr_pos_peak_ptile(sig_dif_msn~= 0);
sig_msn_pos_lfr_peak_ptile = msn_near_lfr_pos_peak_ptile(sig_dif_msn~= 0);
sig_msn_pos_hfr_peak_ratio = msn_near_hfr_pos_peak_ratio(sig_dif_msn~= 0);
sig_msn_pos_lfr_peak_ratio = msn_near_lfr_pos_peak_ratio(sig_dif_msn~= 0);
sig_msn_neg_hfr_peak_ptile = msn_near_hfr_neg_peak_ptile(sig_dif_msn~= 0);
sig_msn_neg_lfr_peak_ptile = msn_near_lfr_neg_peak_ptile(sig_dif_msn~= 0);
sig_msn_neg_hfr_peak_ratio = msn_near_hfr_neg_peak_ratio(sig_dif_msn~= 0);
sig_msn_neg_lfr_peak_ratio = msn_near_lfr_neg_peak_ratio(sig_dif_msn~= 0);


%% Play around with plotting
scatter(msn_near_lfr_pos_peak_ptile(sig_dif_msn ==2), msn_near_hfr_pos_peak_ptile(sig_dif_msn ==2))
figure
scatter(msn_near_lfr_neg_peak_ptile(sig_dif_msn ==1), msn_near_hfr_neg_peak_ptile(sig_dif_msn ==1))


%% Plot PPC Value Histograms
fig_5 = figure;
ax1 = subplot(2,1,1);
hold on;
h1 = histogram(msn_near_ppc_maxval,0:0.01:0.2, 'FaceColor','blue','FaceAlpha',0.4);
h2 = histogram(msn_near_lfr_ppc_maxval,0:0.01:0.2, 'FaceColor','red','FaceAlpha',0.4);
h3 = histogram(msn_near_hfr_ppc_maxval,0:0.01:0.2, 'FaceColor','green','FaceAlpha',0.4);
ax1.Title.String = 'MSN';
ax1.XLabel.String = 'Max PPC Value';

ax2 = subplot(2,1,2);
hold on;
h4 = histogram(fsi_near_ppc_maxval,0:0.01:0.2, 'FaceColor','blue','FaceAlpha',0.4);
h5 = histogram(fsi_near_lfr_ppc_maxval,0:0.01:0.2, 'FaceColor','red','FaceAlpha',0.4);
h6 = histogram(fsi_near_hfr_ppc_maxval,0:0.01:0.2, 'FaceColor','green','FaceAlpha',0.4);
ax2.Title.String = 'FSI';
ax2.XLabel.String = 'Max PPC Value';

%% Plot Only Significant PPC Value Histograms
fig_6 = figure;
ax1 = subplot(2,1,1);
hold on;
h1 = histogram(msn_near_ppc_maxval(sig_dif_msn~=0),0:0.01:0.2, 'FaceColor','blue','FaceAlpha',0.4);
h2 = histogram(msn_near_lfr_ppc_maxval(sig_dif_msn~=0),0:0.01:0.2, 'FaceColor','red','FaceAlpha',0.4);
h3 = histogram(msn_near_hfr_ppc_maxval(sig_dif_msn~=0),0:0.01:0.2, 'FaceColor','green','FaceAlpha',0.4);
ax1.Title.String = 'MSN';
ax1.XLabel.String = 'Max PPC Value';

ax2 = subplot(2,1,2);
hold on;
h4 = histogram(fsi_near_ppc_maxval(sig_dif_fsi~=0),0:0.01:0.2, 'FaceColor','blue','FaceAlpha',0.4);
h5 = histogram(fsi_near_lfr_ppc_maxval(sig_dif_fsi~=0),0:0.01:0.2, 'FaceColor','red','FaceAlpha',0.4);
h6 = histogram(fsi_near_hfr_ppc_maxval(sig_dif_fsi~=0),0:0.01:0.2, 'FaceColor','green','FaceAlpha',0.4);
ax2.Title.String = 'FSI';
ax2.XLabel.String = 'Max PPC Value';

%% Plot distribution of Firing Rates
fig_7 = figure;
ax1 = subplot(2,1,1);
hold on
h1 = histogram(fsi_near_hfr_frs - fsi_near_lfr_frs,0:1:10, 'FaceColor','red','FaceAlpha',0.4);
h2 = histogram(msn_near_hfr_frs - msn_near_lfr_frs,0:1:10, 'FaceColor','green','FaceAlpha',0.4);
leg1 = legend({'FSI', 'MSN'});
ax1.Title.String = 'Distribution of Firing Rate Differences';
ax1.XLabel.String = 'HFR - LFR';

ax2 = subplot(2,1,2);
hold on
h3 = histogram(fsi_near_hfr_frs(sig_dif_fsi~=0) - fsi_near_lfr_frs(sig_dif_fsi~=0),0:1:10, 'FaceColor','red','FaceAlpha',0.4);
h4 = histogram(msn_near_hfr_frs(sig_dif_msn~=0) - msn_near_lfr_frs(sig_dif_msn~=0),0:1:10, 'FaceColor','green','FaceAlpha',0.4);
leg2 = legend({'FSI', 'MSN'});
ax2.Title.String = 'Distribution of Firing Rate Differences for Significant cases';
ax2.XLabel.String = 'HFR - LFR';