cd('D:\RandomVstrAnalysis\ft_results');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
lg = [30,65];
hg = [65,100];
pk_thresh = -1;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        
        fsi_labels  = od.label(od.cell_type == 2);
        for iC = 1:length(fsi_labels)
            o_prefix = extractBefore(fsi_labels{iC},'.t');
            
            % Plot On Track stuff first
            fig = figure('WindowState', 'maximized');
            subplot(3,3,1)
            plot(od.fsi_res.onTrack_spec{iC}.sta_time, od.fsi_res.onTrack_spec{iC}.sta_vals);
            title(sprintf('On Track STA, %d spikes', od.fsi_res.onTrack_spec{iC}.spk_count));
            
            if ~od.fsi_res.onTrack_spec{iC}.flag_nansts
                subplot(3,3,2)
                plot(od.fsi_res.onTrack_spec{iC}.freqs, od.fsi_res.onTrack_spec{iC}.sts_vals);
                xlabel('Freqs')
                title('On Track STS');
            end
            
            if ~od.fsi_res.onTrack_spec{iC}.flag_nanppc
                subplot(3,3,3)
                plot(od.fsi_res.onTrack_spec{iC}.freqs, od.fsi_res.onTrack_spec{iC}.ppc);
                xlabel('Freqs')
                title('On Track PPC');
            end
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if ~od.fsi_res.near_spec{iC}.flag_tooFewSpikes && ...
               ~od.fsi_res.near_spec{iC}.flag_nansts && ...
               ~od.fsi_res.near_spec{iC}.flag_nanppc && ...
               ~od.fsi_res.near_spec{iC}.flag_no_subsampling && ...
               ~od.fsi_res.near_lfr_spec{iC}.flag_tooFewSpikes && ...
               ~od.fsi_res.near_lfr_spec{iC}.flag_nansts && ...
               ~od.fsi_res.near_lfr_spec{iC}.flag_nanppc && ...
               ~od.fsi_res.near_lfr_spec{iC}.flag_no_subsampling && ...
               ~od.fsi_res.near_hfr_spec{iC}.flag_tooFewSpikes && ...
               ~od.fsi_res.near_hfr_spec{iC}.flag_nansts && ...
               ~od.fsi_res.near_hfr_spec{iC}.flag_nanppc && ...
               ~od.fsi_res.near_hfr_spec{iC}.flag_no_subsampling && ...
               ~od.fsi_res.near_p1_spec{iC}.flag_tooFewSpikes && ...
               ~od.fsi_res.near_p1_spec{iC}.flag_nansts && ...
               ~od.fsi_res.near_p1_spec{iC}.flag_nanppc && ...
               ~od.fsi_res.near_p1_spec{iC}.flag_no_subsampling && ...
               ~od.fsi_res.near_p2_spec{iC}.flag_tooFewSpikes && ...
               ~od.fsi_res.near_p2_spec{iC}.flag_nansts && ...
               ~od.fsi_res.near_p2_spec{iC}.flag_nanppc && ...
               ~od.fsi_res.near_p2_spec{iC}.flag_no_subsampling
           
           
                near_subsample_ppc_corr = corrcoef(od.fsi_res.near_spec{iC}.ppc ,od.fsi_res.near_spec{iC}.subsampled_ppc);
                near_subsample_ppc_corr = near_subsample_ppc_corr(1,2);
                near_lfr_subsample_ppc_corr = corrcoef(od.fsi_res.near_lfr_spec{iC}.ppc ,od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                near_lfr_subsample_ppc_corr = near_lfr_subsample_ppc_corr(1,2);
                near_hfr_subsample_ppc_corr = corrcoef(od.fsi_res.near_hfr_spec{iC}.ppc ,od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                near_hfr_subsample_ppc_corr = near_hfr_subsample_ppc_corr(1,2);
                near_subsample_sts_corr = corrcoef(od.fsi_res.near_spec{iC}.sts_vals ,od.fsi_res.near_spec{iC}.subsampled_sts);
                near_subsample_sts_corr = near_subsample_sts_corr(1,2);
                near_lfr_subsample_sts_corr = corrcoef(od.fsi_res.near_lfr_spec{iC}.sts_vals ,od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                near_lfr_subsample_sts_corr = near_lfr_subsample_sts_corr(1,2);
                near_hfr_subsample_sts_corr = corrcoef(od.fsi_res.near_hfr_spec{iC}.sts_vals ,od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                near_hfr_subsample_sts_corr = near_hfr_subsample_sts_corr(1,2);
                
                % Skip if any subsampled correlations are less than 0.9
                if near_subsample_ppc_corr < 0.9 || near_lfr_subsample_ppc_corr < 0.9 || near_hfr_subsample_ppc_corr < 0.9 || ...
                    near_subsample_sts_corr < 0.9 || near_lfr_subsample_sts_corr < 0.9 || near_hfr_subsample_sts_corr < 0.9
                    close all;
                    continue
                end
                
                % Plot STA
                subplot(3,3,4)
                plot(od.fsi_res.near_spec{iC}.sta_time, od.fsi_res.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                flag_leg = true;
                this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.fsi_res.near_spec{iC}.spk_count);
                plot(od.fsi_res.near_spec{iC}.sta_time, od.fsi_res.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.fsi_res.near_lfr_spec{iC}.spk_count);
                plot(od.fsi_res.near_spec{iC}.sta_time, od.fsi_res.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.fsi_res.near_hfr_spec{iC}.spk_count);
                legend(this_legend, 'Location', 'northwest')
                title('Near Reward STA');
                
                flag_near_lg_sts_peak = true;
                flag_near_hg_sts_peak = true;
                flag_near_lg_ppc_peak = true;
                flag_near_hg_ppc_peak = true;
                flag_near_lfr_lg_sts_peak = true;
                flag_near_lfr_hg_sts_peak = true;
                flag_near_lfr_lg_ppc_peak = true;
                flag_near_lfr_hg_ppc_peak = true;
                flag_near_hfr_lg_sts_peak = true;
                flag_near_hfr_hg_sts_peak = true;
                flag_near_hfr_lg_ppc_peak = true;
                flag_near_hfr_hg_ppc_peak = true;
                
                % Peaks in near_spec
                % Find Low Gamma STS Peak
                lf = find(od.fsi_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_lg_sts_peak = false;
                    end
                else
                    flag_near_lg_sts_peak = false; 
                end  
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lg_ppc_peak = false;
                    end
                else
                    flag_near_lg_ppc_peak = false; 
                end
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_hg_sts_peak = false;
                    end
                else
                    flag_near_hg_sts_peak = false; 
                end
                % Find High Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hg_ppc_peak = false;
                    end
                else
                    flag_near_hg_ppc_peak = false; 
                end
                
                % Peaks in near_lfr_spec
                % Find Low Gamma STS Peak
                lf = find(od.fsi_res.near_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_lfr_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_lg_sts_peak = false;
                    end
                else
                    flag_near_lfr_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_lg_ppc_peak = false;
                    end
                else
                    flag_near_lfr_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_lfr_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_hg_sts_peak = false;
                    end
                else
                    flag_near_lfr_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_hg_ppc_peak = false;
                    end
                else
                    flag_near_lfr_hg_ppc_peak = false; 
                end
                
                % Peaks in near_hfr_spec
                % Find Low Gamma STS Peak
                lf = find(od.fsi_res.near_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_hfr_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_lg_sts_peak = false;
                    end
                else
                    flag_near_hfr_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_lg_ppc_peak = false;
                    end
                else
                    flag_near_hfr_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_hfr_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_hfr_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_hg_sts_peak = false;
                    end
                else
                    flag_near_hfr_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_hg_ppc_peak = false;
                    end
                else
                    flag_near_hfr_hg_ppc_peak = false; 
                end
                
                all_sts_peaks_in_lg = flag_near_lg_sts_peak && ...
                    flag_near_hfr_lg_sts_peak && ...
                    flag_near_lfr_lg_sts_peak;                
                all_sts_peaks_in_hg = flag_near_hg_sts_peak && ...
                    flag_near_hfr_hg_sts_peak && ...
                    flag_near_lfr_hg_sts_peak;
                all_ppc_peaks_in_lg = flag_near_lg_ppc_peak && ...
                    flag_near_hfr_lg_ppc_peak && ...
                    flag_near_lfr_lg_ppc_peak;                
                all_ppc_peaks_in_hg = flag_near_hg_ppc_peak && ...
                    flag_near_hfr_hg_ppc_peak && ...
                    flag_near_lfr_hg_ppc_peak;

                % Skip plotting only if no peaks in either hg or lg
                
                if ~(all_sts_peaks_in_lg && all_ppc_peaks_in_lg) && ...
                      ~(all_sts_peaks_in_hg && all_ppc_peaks_in_hg)  
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
                near_lfr_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                near_lfr_sts_corr = near_lfr_sts_corr(1,2);
                near_hfr_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                near_hfr_sts_corr = near_hfr_sts_corr(1,2);
                near_lfr_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                near_lfr_ppc_corr = near_lfr_ppc_corr(1,2);
                near_hfr_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                near_hfr_ppc_corr = near_hfr_ppc_corr(1,2);
              
                % if peaks exist in both low gamma and high gamma range
                if all_sts_peaks_in_lg && all_ppc_peaks_in_lg && ...
                    all_sts_peaks_in_hg && all_ppc_peaks_in_hg ...
                    
                    % Plot STS
                    subplot(3,3,5)
                    this_legend = {};
                    hold on;
                    plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_sts, 'blue');
                    near_lg_sts_pk = od.fsi_res.near_spec{iC}.freqs(lg(1)+near_lg_sts_pk-1);
                    q0 = xline(near_lg_sts_pk, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hg_sts_pk = od.fsi_res.near_spec{iC}.freqs(hg(1)+near_hg_sts_pk-1);
                    q0 = xline(near_hg_sts_pk, '--blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lg_sts_pk, near_hg_sts_pk);
                    plot(od.fsi_res.near_lfr_spec{iC}.freqs, od.fsi_res.near_lfr_spec{iC}.subsampled_sts, 'red');
                    near_lfr_lg_sts_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_sts_pk-1);
                    q0 = xline(near_lfr_lg_sts_pk, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_lfr_hg_sts_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_sts_pk-1);
                    q0 = xline(near_lfr_hg_sts_pk, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lfr_lg_sts_pk, near_lfr_hg_sts_pk);
                    plot(od.fsi_res.near_hfr_spec{iC}.freqs, od.fsi_res.near_hfr_spec{iC}.subsampled_sts, 'green');
                    near_hfr_lg_sts_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_sts_pk-1);
                    q0 = xline(near_hfr_lg_sts_pk, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hfr_hg_sts_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_sts_pk-1);
                    q0 = xline(near_hfr_hg_sts_pk, '--green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('HFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_hfr_lg_sts_pk, near_hfr_hg_sts_pk);
                    sts_text = this_legend;  
%                     legend(this_legend, 'Location', 'northwest')
                    xlabel('Freqs')
                    title('Near Reward STS');
                    
                    
                    % Plot PPC
                    subplot(3,3,6)
                    this_legend = {};
                    hold on;
                    plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_ppc, 'blue');
                    near_lg_ppc_pk = od.fsi_res.near_spec{iC}.freqs(lg(1)+near_lg_ppc_pk-1);
                    q0 = xline(near_lg_ppc_pk, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hg_ppc_pk = od.fsi_res.near_spec{iC}.freqs(hg(1)+near_hg_ppc_pk-1);
                    q0 = xline(near_hg_ppc_pk, '--blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lg_ppc_pk, near_hg_ppc_pk);
                    plot(od.fsi_res.near_lfr_spec{iC}.freqs, od.fsi_res.near_lfr_spec{iC}.subsampled_ppc, 'red');
                    near_lfr_lg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
                    q0 = xline(near_lfr_lg_ppc_pk, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_lfr_hg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
                    q0 = xline(near_lfr_hg_ppc_pk, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lfr_lg_ppc_pk, near_lfr_hg_ppc_pk);
                    plot(od.fsi_res.near_hfr_spec{iC}.freqs, od.fsi_res.near_hfr_spec{iC}.subsampled_ppc, 'green');
                    near_hfr_lg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
                    q0 = xline(near_hfr_lg_ppc_pk, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hfr_hg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);
                    q0 = xline(near_hfr_hg_ppc_pk, '--green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('HFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_hfr_lg_ppc_pk, near_hfr_hg_ppc_pk);
                    ppc_text = this_legend;
%                     legend(this_legend, 'Location', 'northwest')
                    xlabel('Freqs')
                    title('Near Reward PPC');
                    
                   % Relevant info as text
                    subplot(3,3,7)
                    near_text = sprintf("Near Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz\n ", ...
                        near_spec_range(1), near_spec_range(2), near_spec_mean);
                    near_text2 = sprintf("MSN Near spike dist: ");
                    near_text3 = sprintf("%d ", od.msn_near_dist(1:end));
                    near_lfr_text = sprintf("LFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_lfr_spec_range(1), near_lfr_spec_range(2), near_lfr_spec_mean);
                    near_lfr_text2 = sprintf("MSN LFR spike dist:");
                    near_lfr_text3 = sprintf("%d ", od.msn_near_lfr_dist(1:end));
                    near_hfr_text = sprintf("HFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_hfr_spec_range(1), near_hfr_spec_range(2), near_hfr_spec_mean); 
                    near_hfr_text2 = sprintf("MSN HFR spike dist: ");
                    near_hfr_text3 = sprintf("%d ", od.msn_near_hfr_dist(1:end));
                    text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.65, near_text2, 'Color', 'blue', 'FontSize', 10);
                    text(0.4, 0.65, near_text3, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.5, near_lfr_text, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.35, near_lfr_text2, 'Color', 'red', 'FontSize', 10);
                    text(0.4, 0.35, near_lfr_text3, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.2, near_hfr_text, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.05, near_hfr_text2, 'Color', 'green', 'FontSize', 10);
                    text(0.4, 0.05, near_hfr_text3, 'Color', 'green', 'FontSize', 10);
                    axis off;
                    
                    subplot(3,3,8)
                    sts_near_text = sprintf('Correlation between unsampled and subsampled near STS: %.2f', near_subsample_sts_corr);
                    sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
                    sts_lfr_text2 = sprintf('Correlation between unsampled and subsampled LFR STS: %.2f', near_lfr_subsample_sts_corr);
                    sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
                    sts_hfr_text2 = sprintf('Correlation between unsampled and subsampled HFR STS: %.2f', near_hfr_subsample_sts_corr);
                    text(0, 0.9, sts_text{1}, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.8, sts_near_text, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.7, sts_text{2}, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.6, sts_lfr_text, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.5, sts_lfr_text2, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.4, sts_text{3}, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.3, sts_hfr_text, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.2, sts_hfr_text2, 'Color', 'green', 'FontSize', 10);
                    axis off;
                    
                    subplot(3,3,9)
                    ppc_near_text = sprintf('Correlation between unsampled and subsampled near PPC: %.2f', near_subsample_ppc_corr);
                    ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
                    ppc_lfr_text2 = sprintf('Correlation between unsampled and subsampled LFR PPC: %.2f', near_lfr_subsample_ppc_corr);
                    ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
                    ppc_hfr_text2 = sprintf('Correlation between unsampled and subsampled HFR PPC: %.2f', near_hfr_subsample_ppc_corr);
                    text(0, 0.9, ppc_text{1}, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.8, ppc_near_text, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.7, ppc_text{2}, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.6, ppc_lfr_text, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.5, ppc_lfr_text2, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.4, ppc_text{3}, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.3, ppc_hfr_text, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.2, ppc_hfr_text2, 'Color', 'green', 'FontSize', 10);
                    axis off;
                    
                    % Indicate lg category of cell in suffix
                    if near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk > 10
                        o_prefix = cat(2, o_prefix, '_lg_ppc_g10');
                    elseif near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk < -10
                        o_prefix = cat(2, o_prefix, '_lg_ppc_gn10');
                    else
                        o_prefix = cat(2, o_prefix, '_lg_ppc_u10');
                    end
                    
                    % Indicate hg category of cell in suffix
                    if near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk > 10
                        o_prefix = cat(2, o_prefix, '_hg_ppc_g10');
                    elseif near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk < -10
                        o_prefix = cat(2, o_prefix, '_hg_ppc_gn10');
                    else
                        o_prefix = cat(2, o_prefix, '_hg_ppc_u10');
                    end
                    o_name = cat(2, o_prefix, '_FSI');
%                 % If peaks exist only in low gamma    
%                 elseif all_sts_peaks_in_lg && all_ppc_peaks_in_lg
%                     % Plot STS
%                     subplot(3,3,5)
%                     this_legend = {};
%                     hold on;
%                     plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_sts, 'blue');
%                     near_lg_sts_pk = od.fsi_res.near_spec{iC}.freqs(lg(1)+near_lg_sts_pk-1);
%                     q0 = xline(near_lg_sts_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials lg peak: %2.2f Hz ',near_lg_sts_pk);
%                     plot(od.fsi_res.near_lfr_spec{iC}.freqs, od.fsi_res.near_lfr_spec{iC}.subsampled_sts, 'red');
%                     near_lfr_lg_sts_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_sts_pk-1);
%                     q0 = xline(near_lfr_lg_sts_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials lg peak: %2.2f Hz', near_lfr_lg_sts_pk);
%                     plot(od.fsi_res.near_hfr_spec{iC}.freqs, od.fsi_res.near_hfr_spec{iC}.subsampled_sts, 'green');
%                     near_hfr_lg_sts_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_sts_pk-1);
%                     q0 = xline(near_hfr_lg_sts_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials lg peak: %2.2f Hz', near_hfr_lg_sts_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     sts_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward STS');
%                     
%                     % Plot PPC
%                     subplot(3,3,6)
%                     this_legend = {};
%                     hold on;
%                     plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_ppc, 'blue');
%                     near_lg_ppc_pk = od.fsi_res.near_spec{iC}.freqs(lg(1)+near_lg_ppc_pk-1);
%                     q0 = xline(near_lg_ppc_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials lg peak: %2.2f Hz', near_lg_ppc_pk);
%                     plot(od.fsi_res.near_lfr_spec{iC}.freqs, od.fsi_res.near_lfr_spec{iC}.subsampled_ppc, 'red');
%                     near_lfr_lg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
%                     q0 = xline(near_lfr_lg_ppc_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials lg peak: %2.2f Hz', near_lfr_lg_ppc_pk);
%                     plot(od.fsi_res.near_hfr_spec{iC}.freqs, od.fsi_res.near_hfr_spec{iC}.subsampled_ppc, 'green');
%                     near_hfr_lg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
%                     q0 = xline(near_hfr_lg_ppc_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials lg peak: %2.2f Hz', near_hfr_lg_ppc_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     ppc_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward PPC');
%                     
%                     % Relevant info as text
%                     subplot(3,3,7)
%                     near_text = sprintf("Near Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_spec_range(1), near_spec_range(2), near_spec_mean);
%                     near_lfr_text = sprintf("LFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_lfr_spec_range(1), near_lfr_spec_range(2), near_lfr_spec_mean);
%                     near_hfr_text = sprintf("HFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_hfr_spec_range(1), near_hfr_spec_range(2), near_hfr_spec_mean);
%                     text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.5, near_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.2, near_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     axis off;
%                     
%                     subplot(3,3,8)
%                     sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
%                     sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
%                     text(0, 0.8, sts_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, sts_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, sts_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, sts_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, sts_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     sts_hg_flag_text = sprintf("flag_near_lfr_hg_sts_peak: %d\nflag_near_hfr_hg_sts_peak: %d ", ...
%                         flag_near_lfr_hg_sts_peak, flag_near_hfr_hg_sts_peak);
%                     text(0, 0.05, sts_hg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
%                     
%                     subplot(3,3,9)
%                     ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
%                     ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
%                     text(0, 0.8, ppc_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, ppc_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, ppc_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, ppc_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, ppc_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     ppc_hg_flag_text = sprintf("flag_near_lfr_hg_ppc_peak: %d\nflag_near_hfr_hg_ppc_peak: %d ", ...
%                         flag_near_lfr_hg_ppc_peak, flag_near_hfr_hg_ppc_peak);
%                     text(0, 0.05, ppc_hg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
%                     
%                     % Indicate lg category of cell in suffix
%                     if near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk > 10
%                         o_prefix = cat(2, o_prefix, '_lg_ppc_g10');
%                     elseif near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk < -10
%                         o_prefix = cat(2, o_prefix, '_lg_ppc_gn10');
%                     else
%                         o_prefix = cat(2, o_prefix, '_lg_ppc_u10');
%                     end
%                     
%                     % Indicate hg category of cell in suffix
%                     o_prefix = cat(2, o_prefix, '_hg_ppc_NA');   
%                     o_name = cat(2, o_prefix, '_FSI');
%                 % if peaks exist only in high gamma
%                 elseif all_sts_peaks_in_hg && all_ppc_peaks_in_hg
%                     % Plot STS
%                     subplot(3,3,5)
%                     this_legend = {};
%                     hold on;
%                     plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_sts, 'blue');
%                     near_hg_sts_pk = od.fsi_res.near_spec{iC}.freqs(hg(1)+near_hg_sts_pk-1);
%                     q0 = xline(near_hg_sts_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials hg peak: %2.2f Hz ',near_hg_sts_pk);
%                     plot(od.fsi_res.near_lfr_spec{iC}.freqs, od.fsi_res.near_lfr_spec{iC}.subsampled_sts, 'red');
%                     near_lfr_hg_sts_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_sts_pk-1);
%                     q0 = xline(near_lfr_hg_sts_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials hg peak: %2.2f Hz', near_lfr_hg_sts_pk);
%                     plot(od.fsi_res.near_hfr_spec{iC}.freqs, od.fsi_res.near_hfr_spec{iC}.subsampled_sts, 'green');
%                     near_hfr_hg_sts_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_sts_pk-1);
%                     q0 = xline(near_hfr_hg_sts_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials hg peak: %2.2f Hz', near_hfr_hg_sts_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     sts_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward STS');
%                     
%                     % Plot PPC
%                     subplot(3,3,6)
%                     this_legend = {};
%                     hold on;
%                     plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_ppc, 'blue');
%                     near_hg_ppc_pk = od.fsi_res.near_spec{iC}.freqs(hg(1)+near_hg_ppc_pk-1);
%                     q0 = xline(near_hg_ppc_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials hg peak: %2.2f Hz', near_hg_ppc_pk);
%                     plot(od.fsi_res.near_lfr_spec{iC}.freqs, od.fsi_res.near_lfr_spec{iC}.subsampled_ppc, 'red');
%                     near_lfr_hg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
%                     q0 = xline(near_lfr_hg_ppc_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials hg peak: %2.2f Hz', near_lfr_hg_ppc_pk);
%                     plot(od.fsi_res.near_hfr_spec{iC}.freqs, od.fsi_res.near_hfr_spec{iC}.subsampled_ppc, 'green');
%                     near_hfr_hg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);
%                     q0 = xline(near_hfr_hg_ppc_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials hg peak: %2.2f Hz', near_hfr_hg_ppc_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     ppc_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward PPC');
%                     
%                      % Relevant info as text
%                     subplot(3,3,7)
%                     near_text = sprintf("Near Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_spec_range(1), near_spec_range(2), near_spec_mean);
%                     near_lfr_text = sprintf("LFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_lfr_spec_range(1), near_lfr_spec_range(2), near_lfr_spec_mean);
%                     near_hfr_text = sprintf("HFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_hfr_spec_range(1), near_hfr_spec_range(2), near_hfr_spec_mean);
%                     text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.5, near_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.2, near_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     axis off;
%                     
%                     subplot(3,3,8)
%                     sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
%                     sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
%                     text(0, 0.8, sts_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, sts_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, sts_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, sts_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, sts_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     sts_lg_flag_text = sprintf("flag_near_lfr_lg_sts_peak: %d\nflag_near_hfr_lg_sts_peak: %d ", ...
%                         flag_near_lfr_lg_sts_peak, flag_near_hfr_lg_sts_peak);
%                     text(0, 0.05, sts_lg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
%                     
%                     subplot(3,3,9)
%                     ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
%                     ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
%                     text(0, 0.8, ppc_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, ppc_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, ppc_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, ppc_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, ppc_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     ppc_lg_flag_text = sprintf("flag_near_lfr_lg_ppc_peak: %d\nflag_near_hfr_lg_ppc_peak: %d ", ...
%                         flag_near_lfr_lg_ppc_peak, flag_near_hfr_lg_ppc_peak);
%                     text(0, 0.05, ppc_lg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
% 
%                     % Indicate hg category of cell in suffix
%                     o_prefix = cat(2, o_prefix, '_lg_ppc_NA');
%                     
%                     % Indicate hg category of cell in suffix
%                     if near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk > 10
%                         o_prefix = cat(2, o_prefix, '_hg_ppc_g10');
%                     elseif near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk < -10
%                         o_prefix = cat(2, o_prefix, '_hg_ppc_gn10');
%                     else
%                         o_prefix = cat(2, o_prefix, '_hg_ppc_u10');
%                     end 
%                     o_name = cat(2, o_prefix, '_FSI');
                end
            else
                close all;
                continue
            end
            WriteFig(fig,o_name,1);
            close all;  
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        for iC = 1:length(msn_labels)
            o_prefix = extractBefore(msn_labels{iC},'.t');
            
            % Plot On Track stuff first
            fig = figure('WindowState', 'maximized');
            subplot(3,3,1)
            plot(od.msn_res.onTrack_spec{iC}.sta_time, od.msn_res.onTrack_spec{iC}.sta_vals);
            title(sprintf('On Track STA, %d spikes', od.msn_res.onTrack_spec{iC}.spk_count));
            
            if ~od.msn_res.onTrack_spec{iC}.flag_nansts
                subplot(3,3,2)
                plot(od.msn_res.onTrack_spec{iC}.freqs, od.msn_res.onTrack_spec{iC}.sts_vals);
                xlabel('Freqs')
                title('On Track STS');
            end
            
            if ~od.msn_res.onTrack_spec{iC}.flag_nanppc
                subplot(3,3,3)
                plot(od.msn_res.onTrack_spec{iC}.freqs, od.msn_res.onTrack_spec{iC}.ppc);
                xlabel('Freqs')
                title('On Track PPC');
            end
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if ~od.msn_res.near_spec{iC}.flag_tooFewSpikes && ...
               ~od.msn_res.near_spec{iC}.flag_nansts && ...
               ~od.msn_res.near_spec{iC}.flag_nanppc && ...
               ~od.msn_res.near_lfr_spec{iC}.flag_tooFewSpikes && ...
               ~od.msn_res.near_lfr_spec{iC}.flag_nansts && ...
               ~od.msn_res.near_lfr_spec{iC}.flag_nanppc && ...
               ~od.msn_res.near_hfr_spec{iC}.flag_tooFewSpikes && ...
               ~od.msn_res.near_hfr_spec{iC}.flag_nansts && ...
               ~od.msn_res.near_hfr_spec{iC}.flag_nanppc && ...
               ~od.msn_res.near_p1_spec{iC}.flag_tooFewSpikes && ...
               ~od.msn_res.near_p1_spec{iC}.flag_nansts && ...
               ~od.msn_res.near_p1_spec{iC}.flag_nanppc && ...
               ~od.msn_res.near_p2_spec{iC}.flag_tooFewSpikes && ...
               ~od.msn_res.near_p2_spec{iC}.flag_nansts && ...
               ~od.msn_res.near_p2_spec{iC}.flag_nanppc
                
                % Plot STA
                subplot(3,3,4)
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                flag_leg = true;
                this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.msn_res.near_spec{iC}.spk_count);
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.msn_res.near_lfr_spec{iC}.spk_count);
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.msn_res.near_hfr_spec{iC}.spk_count);
                legend(this_legend, 'Location', 'northwest')
                title('Near Reward STA');
                
                flag_near_lg_sts_peak = true;
                flag_near_hg_sts_peak = true;
                flag_near_lg_ppc_peak = true;
                flag_near_hg_ppc_peak = true;
                flag_near_lfr_lg_sts_peak = true;
                flag_near_lfr_hg_sts_peak = true;
                flag_near_lfr_lg_ppc_peak = true;
                flag_near_lfr_hg_ppc_peak = true;
                flag_near_hfr_lg_sts_peak = true;
                flag_near_hfr_hg_sts_peak = true;
                flag_near_hfr_lg_ppc_peak = true;
                flag_near_hfr_hg_ppc_peak = true;
                
                % Peaks in near_spec
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_lg_sts_peak = false;
                    end
                else
                    flag_near_lg_sts_peak = false; 
                end  
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lg_ppc_peak = false;
                    end
                else
                    flag_near_lg_ppc_peak = false; 
                end
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_hg_sts_peak = false;
                    end
                else
                    flag_near_hg_sts_peak = false; 
                end
                % Find High Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hg_ppc_peak = false;
                    end
                else
                    flag_near_hg_ppc_peak = false; 
                end
                
                % Peaks in near_lfr_spec
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_lfr_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_lfr_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_lg_sts_peak = false;
                    end
                else
                    flag_near_lfr_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_lfr_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_lfr_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_lg_ppc_peak = false;
                    end
                else
                    flag_near_lfr_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_lfr_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_lfr_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_hg_sts_peak = false;
                    end
                else
                    flag_near_lfr_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_lfr_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_lfr_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_lfr_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_lfr_hg_ppc_peak = false;
                    end
                else
                    flag_near_lfr_hg_ppc_peak = false; 
                end
                
                % Peaks in near_hfr_spec
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_hfr_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_hfr_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_lg_sts_peak = false;
                    end
                else
                    flag_near_hfr_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_hfr_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_hfr_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_lg_ppc_peak = false;
                    end
                else
                    flag_near_hfr_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_hfr_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_hfr_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_hfr_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_hg_sts_peak = false;
                    end
                else
                    flag_near_hfr_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_hfr_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_hfr_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_hfr_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hfr_hg_ppc_peak = false;
                    end
                else
                    flag_near_hfr_hg_ppc_peak = false; 
                end
                
                all_sts_peaks_in_lg = flag_near_lg_sts_peak && ...
                    flag_near_hfr_lg_sts_peak && ...
                    flag_near_lfr_lg_sts_peak;                
                all_sts_peaks_in_hg = flag_near_hg_sts_peak && ...
                    flag_near_hfr_hg_sts_peak && ...
                    flag_near_lfr_hg_sts_peak;
                all_ppc_peaks_in_lg = flag_near_lg_ppc_peak && ...
                    flag_near_hfr_lg_ppc_peak && ...
                    flag_near_lfr_lg_ppc_peak;                
                all_ppc_peaks_in_hg = flag_near_hg_ppc_peak && ...
                    flag_near_hfr_hg_ppc_peak && ...
                    flag_near_lfr_hg_ppc_peak;

                % Skip plotting only if no peaks in either hg or lg
                
                if ~(all_sts_peaks_in_lg && all_ppc_peaks_in_lg) && ...
                      ~(all_sts_peaks_in_hg && all_ppc_peaks_in_hg)  
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
                near_spec_mean = sum(od.msn_res.near_spec{iC}.trial_spk_count(nz_trials))/sum(od.msn_res.near_spec{iC}.trial_spk_count(nz_trials)'./od.msn_res.near_spec{iC}.mfr(nz_trials));
                near_spec_range = [min(od.msn_res.near_spec{iC}.mfr(nz_trials)), max(od.msn_res.near_spec{iC}.mfr(nz_trials))];
                near_lfr_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,od.msn_res.near_lfr_spec{iC}.sts_vals);
                near_lfr_sts_corr = near_lfr_sts_corr(1,2);
                near_hfr_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,od.msn_res.near_hfr_spec{iC}.sts_vals);
                near_hfr_sts_corr = near_hfr_sts_corr(1,2);
                near_lfr_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,od.msn_res.near_lfr_spec{iC}.ppc);
                near_lfr_ppc_corr = near_lfr_ppc_corr(1,2);
                near_hfr_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,od.msn_res.near_hfr_spec{iC}.ppc);
                near_hfr_ppc_corr = near_hfr_ppc_corr(1,2);
                
                % if peaks exist in both low gamma and high gamma range
                if all_sts_peaks_in_lg && all_ppc_peaks_in_lg && ...
                    all_sts_peaks_in_hg && all_ppc_peaks_in_hg ...
                    
                    % Plot STS
                    subplot(3,3,5)
                    this_legend = {};
                    hold on;
                    plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.sts_vals, 'blue');
                    near_lg_sts_pk = od.msn_res.near_spec{iC}.freqs(lg(1)+near_lg_sts_pk-1);
                    q0 = xline(near_lg_sts_pk, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hg_sts_pk = od.msn_res.near_spec{iC}.freqs(hg(1)+near_hg_sts_pk-1);
                    q0 = xline(near_hg_sts_pk, '--blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lg_sts_pk, near_hg_sts_pk);
                    plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.sts_vals, 'red');
                    near_lfr_lg_sts_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_sts_pk-1);
                    q0 = xline(near_lfr_lg_sts_pk, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_lfr_hg_sts_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_sts_pk-1);
                    q0 = xline(near_lfr_hg_sts_pk, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lfr_lg_sts_pk, near_lfr_hg_sts_pk);
                    plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.sts_vals, 'green');
                    near_hfr_lg_sts_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_sts_pk-1);
                    q0 = xline(near_hfr_lg_sts_pk, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hfr_hg_sts_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_sts_pk-1);
                    q0 = xline(near_hfr_hg_sts_pk, '--green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('HFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_hfr_lg_sts_pk, near_hfr_hg_sts_pk);
                    sts_text = this_legend;  
%                     legend(this_legend, 'Location', 'northwest')
                    xlabel('Freqs')
                    title('Near Reward STS');
                    
                    
                    % Plot PPC
                    subplot(3,3,6)
                    this_legend = {};
                    hold on;
                    plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.ppc, 'blue');
                    near_lg_ppc_pk = od.msn_res.near_spec{iC}.freqs(lg(1)+near_lg_ppc_pk-1);
                    q0 = xline(near_lg_ppc_pk, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hg_ppc_pk = od.msn_res.near_spec{iC}.freqs(hg(1)+near_hg_ppc_pk-1);
                    q0 = xline(near_hg_ppc_pk, '--blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lg_ppc_pk, near_hg_ppc_pk);
                    plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.ppc, 'red');
                    near_lfr_lg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
                    q0 = xline(near_lfr_lg_ppc_pk, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_lfr_hg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
                    q0 = xline(near_lfr_hg_ppc_pk, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_lfr_lg_ppc_pk, near_lfr_hg_ppc_pk);
                    plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.ppc, 'green');
                    near_hfr_lg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
                    q0 = xline(near_hfr_lg_ppc_pk, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    near_hfr_hg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);
                    q0 = xline(near_hfr_hg_ppc_pk, '--green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                       ('HFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                       near_hfr_lg_ppc_pk, near_hfr_hg_ppc_pk);
                    ppc_text = this_legend;
%                     legend(this_legend, 'Location', 'northwest')
                    xlabel('Freqs')
                    title('Near Reward PPC');
                    
                    % Relevant info as text
                    subplot(3,3,7)
                    near_text = sprintf("Near Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_spec_range(1), near_spec_range(2), near_spec_mean);
                    near_lfr_text = sprintf("LFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_lfr_spec_range(1), near_lfr_spec_range(2), near_lfr_spec_mean);
                    near_hfr_text = sprintf("HFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_hfr_spec_range(1), near_hfr_spec_range(2), near_hfr_spec_mean);
                    text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.5, near_lfr_text, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.2, near_hfr_text, 'Color', 'green', 'FontSize', 10);
                    axis off;
                    
                    subplot(3,3,8)
                    sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
                    sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
                    text(0, 0.8, sts_text{1}, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.65, sts_text{2}, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.5, sts_lfr_text, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.35, sts_text{3}, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.2, sts_hfr_text, 'Color', 'green', 'FontSize', 10);
                    axis off;
                    
                    subplot(3,3,9)
                    ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
                    ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
                    text(0, 0.8, ppc_text{1}, 'Color', 'blue', 'FontSize', 10);
                    text(0, 0.65, ppc_text{2}, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.5, ppc_lfr_text, 'Color', 'red', 'FontSize', 10);
                    text(0, 0.35, ppc_text{3}, 'Color', 'green', 'FontSize', 10);
                    text(0, 0.2, ppc_hfr_text, 'Color', 'green', 'FontSize', 10);
                    axis off;
                    
                    % Indicate lg category of cell in suffix
                    if near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk > 10
                        o_prefix = cat(2, o_prefix, '_lg_ppc_g10');
                    elseif near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk < -10
                        o_prefix = cat(2, o_prefix, '_lg_ppc_gn10');
                    else
                        o_prefix = cat(2, o_prefix, '_lg_ppc_u10');
                    end
                    
                    % Indicate hg category of cell in suffix
                    if near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk > 10
                        o_prefix = cat(2, o_prefix, '_hg_ppc_g10');
                    elseif near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk < -10
                        o_prefix = cat(2, o_prefix, '_hg_ppc_gn10');
                    else
                        o_prefix = cat(2, o_prefix, '_hg_ppc_u10');
                    end
                    o_name = cat(2, o_prefix, '_MSN');
%                 % If peaks exist only in low gamma    
%                 elseif all_sts_peaks_in_lg && all_ppc_peaks_in_lg
%                     % Plot STS
%                     subplot(3,3,5)
%                     this_legend = {};
%                     hold on;
%                     plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.sts_vals, 'blue');
%                     near_lg_sts_pk = od.msn_res.near_spec{iC}.freqs(lg(1)+near_lg_sts_pk-1);
%                     q0 = xline(near_lg_sts_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials lg peak: %2.2f Hz ',near_lg_sts_pk);
%                     plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.sts_vals, 'red');
%                     near_lfr_lg_sts_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_sts_pk-1);
%                     q0 = xline(near_lfr_lg_sts_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials lg peak: %2.2f Hz', near_lfr_lg_sts_pk);
%                     plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.sts_vals, 'green');
%                     near_hfr_lg_sts_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_sts_pk-1);
%                     q0 = xline(near_hfr_lg_sts_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials lg peak: %2.2f Hz', near_hfr_lg_sts_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     sts_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward STS');
%                     
%                     % Plot PPC
%                     subplot(3,3,6)
%                     this_legend = {};
%                     hold on;
%                     plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.ppc, 'blue');
%                     near_lg_ppc_pk = od.msn_res.near_spec{iC}.freqs(lg(1)+near_lg_ppc_pk-1);
%                     q0 = xline(near_lg_ppc_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials lg peak: %2.2f Hz', near_lg_ppc_pk);
%                     plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.ppc, 'red');
%                     near_lfr_lg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
%                     q0 = xline(near_lfr_lg_ppc_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials lg peak: %2.2f Hz', near_lfr_lg_ppc_pk);
%                     plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.ppc, 'green');
%                     near_hfr_lg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
%                     q0 = xline(near_hfr_lg_ppc_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials lg peak: %2.2f Hz', near_hfr_lg_ppc_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     ppc_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward PPC');
%                     
%                     % Relevant info as text
%                     subplot(3,3,7)
%                     near_text = sprintf("Near Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_spec_range(1), near_spec_range(2), near_spec_mean);
%                     near_lfr_text = sprintf("LFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_lfr_spec_range(1), near_lfr_spec_range(2), near_lfr_spec_mean);
%                     near_hfr_text = sprintf("HFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_hfr_spec_range(1), near_hfr_spec_range(2), near_hfr_spec_mean);
%                     text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.5, near_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.2, near_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     axis off;
%                     
%                     subplot(3,3,8)
%                     sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
%                     sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
%                     text(0, 0.8, sts_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, sts_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, sts_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, sts_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, sts_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     sts_hg_flag_text = sprintf("flag_near_lfr_hg_sts_peak: %d\nflag_near_hfr_hg_sts_peak: %d ", ...
%                         flag_near_lfr_hg_sts_peak, flag_near_hfr_hg_sts_peak);
%                     text(0, 0.05, sts_hg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
%                     
%                     subplot(3,3,9)
%                     ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
%                     ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
%                     text(0, 0.8, ppc_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, ppc_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, ppc_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, ppc_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, ppc_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     ppc_hg_flag_text = sprintf("flag_near_lfr_hg_ppc_peak: %d\nflag_near_hfr_hg_ppc_peak: %d ", ...
%                         flag_near_lfr_hg_ppc_peak, flag_near_hfr_hg_ppc_peak);
%                     text(0, 0.05, ppc_hg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
%                     
%                     % Indicate lg category of cell in suffix
%                     if near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk > 10
%                         o_prefix = cat(2, o_prefix, '_lg_ppc_g10');
%                     elseif near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk < -10
%                         o_prefix = cat(2, o_prefix, '_lg_ppc_gn10');
%                     else
%                         o_prefix = cat(2, o_prefix, '_lg_ppc_u10');
%                     end
%                     
%                     % Indicate hg category of cell in suffix
%                     o_prefix = cat(2, o_prefix, '_hg_ppc_NA');   
%                     o_name = cat(2, o_prefix, '_MSN');
%                 % if peaks exist only in high gamma
%                 elseif all_sts_peaks_in_hg && all_ppc_peaks_in_hg
%                     % Plot STS
%                     subplot(3,3,5)
%                     this_legend = {};
%                     hold on;
%                     plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.sts_vals, 'blue');
%                     near_hg_sts_pk = od.msn_res.near_spec{iC}.freqs(hg(1)+near_hg_sts_pk-1);
%                     q0 = xline(near_hg_sts_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials hg peak: %2.2f Hz ',near_hg_sts_pk);
%                     plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.sts_vals, 'red');
%                     near_lfr_hg_sts_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_sts_pk-1);
%                     q0 = xline(near_lfr_hg_sts_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials hg peak: %2.2f Hz', near_lfr_hg_sts_pk);
%                     plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.sts_vals, 'green');
%                     near_hfr_hg_sts_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_sts_pk-1);
%                     q0 = xline(near_hfr_hg_sts_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials hg peak: %2.2f Hz', near_hfr_hg_sts_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     sts_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward STS');
%                     
%                     % Plot PPC
%                     subplot(3,3,6)
%                     this_legend = {};
%                     hold on;
%                     plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.ppc, 'blue');
%                     near_hg_ppc_pk = od.msn_res.near_spec{iC}.freqs(hg(1)+near_hg_ppc_pk-1);
%                     q0 = xline(near_hg_ppc_pk, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('All Trials hg peak: %2.2f Hz', near_hg_ppc_pk);
%                     plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.ppc, 'red');
%                     near_lfr_hg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
%                     q0 = xline(near_lfr_hg_ppc_pk, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('LFR Trials hg peak: %2.2f Hz', near_lfr_hg_ppc_pk);
%                     plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.ppc, 'green');
%                     near_hfr_hg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);
%                     q0 = xline(near_hfr_hg_ppc_pk, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                        ('HFR Trials hg peak: %2.2f Hz', near_hfr_hg_ppc_pk);
% %                     legend(this_legend, 'Location', 'northwest')
%                     ppc_text = this_legend;
%                     xlabel('Freqs')
%                     title('Near Reward PPC');
%                     
%                      % Relevant info as text
%                     subplot(3,3,7)
%                     near_text = sprintf("Near Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_spec_range(1), near_spec_range(2), near_spec_mean);
%                     near_lfr_text = sprintf("LFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_lfr_spec_range(1), near_lfr_spec_range(2), near_lfr_spec_mean);
%                     near_hfr_text = sprintf("HFR Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
%                         near_hfr_spec_range(1), near_hfr_spec_range(2), near_hfr_spec_mean);
%                     text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.5, near_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.2, near_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     axis off;
%                     
%                     subplot(3,3,8)
%                     sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
%                     sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
%                     text(0, 0.8, sts_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, sts_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, sts_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, sts_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, sts_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     sts_lg_flag_text = sprintf("flag_near_lfr_lg_sts_peak: %d\nflag_near_hfr_lg_sts_peak: %d ", ...
%                         flag_near_lfr_lg_sts_peak, flag_near_hfr_lg_sts_peak);
%                     text(0, 0.05, sts_lg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
%                     
%                     subplot(3,3,9)
%                     ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
%                     ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
%                     text(0, 0.8, ppc_text{1}, 'Color', 'blue', 'FontSize', 10);
%                     text(0, 0.65, ppc_text{2}, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.5, ppc_lfr_text, 'Color', 'red', 'FontSize', 10);
%                     text(0, 0.35, ppc_text{3}, 'Color', 'green', 'FontSize', 10);
%                     text(0, 0.2, ppc_hfr_text, 'Color', 'green', 'FontSize', 10);
%                     ppc_lg_flag_text = sprintf("flag_near_lfr_lg_ppc_peak: %d\nflag_near_hfr_lg_ppc_peak: %d ", ...
%                         flag_near_lfr_lg_ppc_peak, flag_near_hfr_lg_ppc_peak);
%                     text(0, 0.05, ppc_lg_flag_text, 'Color', 'black', 'FontSize', 10, 'Interpreter', 'none');
%                     axis off;
% 
%                     % Indicate hg category of cell in suffix
%                     o_prefix = cat(2, o_prefix, '_lg_ppc_NA');
%                     
%                     % Indicate hg category of cell in suffix
%                     if near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk > 10
%                         o_prefix = cat(2, o_prefix, '_hg_ppc_g10');
%                     elseif near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk < -10
%                         o_prefix = cat(2, o_prefix, '_hg_ppc_gn10');
%                     else
%                         o_prefix = cat(2, o_prefix, '_hg_ppc_u10');
%                     end 
%                     o_name = cat(2, o_prefix, '_MSN');
                end
            else
                close all;
                continue
            end
            WriteFig(fig,o_name,1);
            close all;  
        end
        
        
        
    end
end