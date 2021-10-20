cd('D:\RandomVstrAnalysis\ft_results');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
lg = [5,30];
hg = [30,100];
pk_thresh = -1;
num_control_splits = 100;

for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        
        fsi_labels  = od.label(od.cell_type == 2);
        for iC = 1:length(fsi_labels)
            o_prefix = extractBefore(fsi_labels{iC},'.t');
            fig = figure('WindowState', 'maximized');
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if isfield(od.fsi_res.near_spec{iC},'flag_no_control_split') && ...
                        ~od.fsi_res.near_spec{iC}.flag_no_control_split          
                near_subsample_ppc_corr = corrcoef(od.fsi_res.near_spec{iC}.ppc ,od.fsi_res.near_spec{iC}.subsampled_ppc);
                near_subsample_ppc_corr = near_subsample_ppc_corr(1,2);
                near_lfr_subsample_ppc_corr = corrcoef(od.fsi_res.near_lfr_spec{iC}.ppc ,od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                near_lfr_subsample_ppc_corr = near_lfr_subsample_ppc_corr(1,2);
                near_hfr_subsample_ppc_corr = corrcoef(od.fsi_res.near_hfr_spec{iC}.ppc ,od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                near_hfr_subsample_ppc_corr = near_hfr_subsample_ppc_corr(1,2);
                near_p1_subsample_ppc_corr = corrcoef(mean(od.fsi_res.near_p1_spec{iC}.ppc,1) ,mean(od.fsi_res.near_p1_spec{iC}.subsampled_ppc,1));
                near_p1_subsample_ppc_corr = near_p1_subsample_ppc_corr(1,2);
                near_p2_subsample_ppc_corr = corrcoef(mean(od.fsi_res.near_p2_spec{iC}.ppc,1) ,mean(od.fsi_res.near_p2_spec{iC}.subsampled_ppc,1));
                near_p2_subsample_ppc_corr = near_p2_subsample_ppc_corr(1,2);
                near_subsample_sts_corr = corrcoef(od.fsi_res.near_spec{iC}.sts_vals ,od.fsi_res.near_spec{iC}.subsampled_sts);
                near_subsample_sts_corr = near_subsample_sts_corr(1,2);
                near_lfr_subsample_sts_corr = corrcoef(od.fsi_res.near_lfr_spec{iC}.sts_vals ,od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                near_lfr_subsample_sts_corr = near_lfr_subsample_sts_corr(1,2);
                near_hfr_subsample_sts_corr = corrcoef(od.fsi_res.near_hfr_spec{iC}.sts_vals ,od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                near_hfr_subsample_sts_corr = near_hfr_subsample_sts_corr(1,2);
                near_p1_subsample_sts_corr = corrcoef(mean(od.fsi_res.near_p1_spec{iC}.sts,1) ,mean(od.fsi_res.near_p1_spec{iC}.subsampled_sts,1));
                near_p1_subsample_sts_corr = near_p1_subsample_sts_corr(1,2);
                near_p2_subsample_sts_corr = corrcoef(mean(od.fsi_res.near_p2_spec{iC}.sts,1) ,mean(od.fsi_res.near_p2_spec{iC}.subsampled_sts,1));
                near_p2_subsample_sts_corr = near_p2_subsample_sts_corr(1,2);
                
%                 % Maybe this step is unnecessary?
%                 % Skip if any subsampled correlations are less than 0.9
%                 if near_subsample_ppc_corr < 0.9 || near_lfr_subsample_ppc_corr < 0.9 || near_hfr_subsample_ppc_corr < 0.9 || ...
%                     near_subsample_sts_corr < 0.9 || near_lfr_subsample_sts_corr < 0.9 || near_hfr_subsample_sts_corr < 0.9 || ...
%                     near_p1_subsample_ppc_corr < 0.9 || near_p2_subsample_ppc_corr < 0.9 || ...
%                     near_p1_subsample_sts_corr < 0.9 || near_p2_subsample_sts_corr < 0.9
%                     close all;
%                     continue
%                 end
                
                % Plot STA
                subplot(3,3,1)
                plot(od.fsi_res.near_spec{iC}.sta_time, od.fsi_res.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                flag_leg = true;
                this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.fsi_res.near_spec{iC}.spk_count);
                plot(od.fsi_res.near_spec{iC}.sta_time, od.fsi_res.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.fsi_res.near_lfr_spec{iC}.spk_count);
                plot(od.fsi_res.near_spec{iC}.sta_time, od.fsi_res.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.fsi_res.near_hfr_spec{iC}.spk_count);
                legend(this_legend, 'Location', 'northwest', 'FontSize', 8)
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
                flag_near_p1_lg_sts_peak = true;
                flag_near_p1_hg_sts_peak = true;
                flag_near_p1_lg_ppc_peak = true;
                flag_near_p1_hg_ppc_peak = true;
                flag_near_p2_lg_sts_peak = true;
                flag_near_p2_hg_sts_peak = true;
                flag_near_p2_lg_ppc_peak = true;
                flag_near_p2_hg_ppc_peak = true;
                near_p1_lg_sts_pk = zeros(1,num_control_splits);
                near_p2_lg_sts_pk = zeros(1,num_control_splits);
                near_p1_hg_sts_pk = zeros(1,num_control_splits);
                near_p2_hg_sts_pk = zeros(1,num_control_splits);
                near_p1_lg_ppc_pk = zeros(1,num_control_splits);
                near_p2_hg_ppc_pk = zeros(1,num_control_splits);
                near_p1_hg_ppc_pk = zeros(1,num_control_splits);
                near_p2_lg_ppc_pk = zeros(1,num_control_splits);
                

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
                
                % Peaks in near_p1_spec
                % Find Low Gamma STS Peak
                lf = find(od.fsi_res.near_p1_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_p1_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_sts(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_lg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_lg_sts_peak = false;
                        end
                    else
                        flag_near_p1_lg_sts_peak = false; 
                    end
                end   
                % Find Low Gamma PPC Peak
                lf = find(od.fsi_res.near_p1_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_p1_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_lg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_lg_ppc_peak = false;
                        end
                    else
                        flag_near_p1_lg_ppc_peak = false; 
                    end
                end
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_p1_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_p1_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_sts(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_hg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_hg_sts_peak = false;
                        end
                    else
                        flag_near_p1_hg_sts_peak = false; 
                    end
                end   
                % Find High Gamma PPC Peak
                lf = find(od.fsi_res.near_p1_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_p1_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_hg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_hg_ppc_peak = false;
                        end
                    else
                        flag_near_p1_hg_ppc_peak = false; 
                    end
                end
                
                % Peaks in near_p2_spec
                % Find Low Gamma STS Peak
                lf = find(od.fsi_res.near_p2_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_p2_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_sts(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_lg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_lg_sts_peak = false;
                        end
                    else
                        flag_near_p2_lg_sts_peak = false; 
                    end
                end   
                % Find Low Gamma PPC Peak
                lf = find(od.fsi_res.near_p2_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_p2_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_lg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_lg_ppc_peak = false;
                        end
                    else
                        flag_near_p2_lg_ppc_peak = false; 
                    end
                end
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_p2_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_p2_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_sts(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_hg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_hg_sts_peak = false;
                        end
                    else
                        flag_near_p2_hg_sts_peak = false; 
                    end
                end   
                % Find High Gamma PPC Peak
                lf = find(od.fsi_res.near_p2_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_p2_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(iP,lf:rf));
                    max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_hg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_hg_ppc_peak = false;
                        end
                    else
                        flag_near_p2_hg_ppc_peak = false; 
                    end
                end
                
                all_sts_peaks_in_lg = flag_near_lg_sts_peak && ...
                    flag_near_hfr_lg_sts_peak && ...
                    flag_near_lfr_lg_sts_peak && ...
                    flag_near_p2_lg_sts_peak && ...
                    flag_near_p1_lg_sts_peak;
                all_sts_peaks_in_hg = flag_near_hg_sts_peak && ...
                    flag_near_hfr_hg_sts_peak && ...
                    flag_near_lfr_hg_sts_peak && ...
                    flag_near_p2_hg_sts_peak && ...
                    flag_near_p1_hg_sts_peak;
                all_ppc_peaks_in_lg = flag_near_lg_ppc_peak && ...
                    flag_near_hfr_lg_ppc_peak && ...
                    flag_near_lfr_lg_ppc_peak && ...
                    flag_near_p2_lg_ppc_peak && ...
                    flag_near_p1_lg_ppc_peak;                              
                all_ppc_peaks_in_hg = flag_near_hg_ppc_peak && ...
                    flag_near_hfr_hg_ppc_peak && ...
                    flag_near_lfr_hg_ppc_peak && ...
                    flag_near_p2_hg_ppc_peak && ...
                    flag_near_p1_hg_ppc_peak;

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
                near_lfr_sts_corr = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts,od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                near_lfr_sts_corr = near_lfr_sts_corr(1,2);
                near_hfr_sts_corr = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts,od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                near_hfr_sts_corr = near_hfr_sts_corr(1,2);
                near_lfr_ppc_corr = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc,od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                near_lfr_ppc_corr = near_lfr_ppc_corr(1,2);
                near_hfr_ppc_corr = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc,od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                near_hfr_ppc_corr = near_hfr_ppc_corr(1,2);
                 
                p1_trials = od.fsi_res.near_spec{iC}.valid_splits;    
                p2_trials = ~p1_trials;
                p1_trials = p1_trials & repmat(nz_trials',100,1);
                p2_trials = p2_trials & repmat(nz_trials',100,1);
                temp_p1_spk = zeros(100,1);
                temp_p1_dur = zeros(100,1);
                temp_p2_spk = zeros(100,1);
                temp_p2_dur = zeros(100,1);
                for iP = 1:100
                    temp_p1_spk(iP) = sum(od.fsi_res.near_spec{iC}.trial_spk_count(p1_trials(iP,:)));
                    temp_p1_dur(iP) = sum(od.fsi_res.near_spec{iC}.trial_spk_count(p1_trials(iP,:))./od.fsi_res.near_spec{iC}.mfr(p1_trials(iP,:))');
                    temp_p2_spk(iP) = sum(od.fsi_res.near_spec{iC}.trial_spk_count(p2_trials(iP,:)));
                    temp_p2_dur(iP) = sum(od.fsi_res.near_spec{iC}.trial_spk_count(p2_trials(iP,:))./od.fsi_res.near_spec{iC}.mfr(p2_trials(iP,:))');
                end
                temp_p1_mfr = temp_p1_spk./temp_p1_dur;
                temp_p2_mfr = temp_p2_spk./temp_p2_dur;
                near_p1_spec_mean = mean(temp_p1_mfr);
                near_p1_spec_range = [min(temp_p1_mfr) max(temp_p1_mfr)];
                near_p2_spec_mean = mean(temp_p2_mfr);
                near_p2_spec_range = [min(temp_p2_mfr) max(temp_p2_mfr)];
                near_p1_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,mean(od.fsi_res.near_p1_spec{iC}.subsampled_sts,1));
                near_p1_sts_corr = near_p1_sts_corr(1,2);
                near_p2_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,mean(od.fsi_res.near_p2_spec{iC}.subsampled_sts,1));
                near_p2_sts_corr = near_p2_sts_corr(1,2);
                near_p1_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,mean(od.fsi_res.near_p1_spec{iC}.subsampled_ppc,1));
                near_p1_ppc_corr = near_p1_ppc_corr(1,2);
                near_p2_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,mean(od.fsi_res.near_p2_spec{iC}.subsampled_ppc,1));
                near_p2_ppc_corr = near_p2_ppc_corr(1,2);
              
                % if peaks exist in both low gamma and high gamma range
                if all_sts_peaks_in_lg && all_ppc_peaks_in_lg && ...
                    all_sts_peaks_in_hg && all_ppc_peaks_in_hg ...
                    

                    % Plot STS
                    subplot(3,3,2)
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
                    x1 = find(od.fsi_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    xlabel('Freqs')
                    title('Near Reward STS');
                   
                    p1_p2_sts = od.fsi_res.near_p1_spec{iC}.subsampled_sts - od.fsi_res.near_p2_spec{iC}.subsampled_sts;
                    p1_p2_ppc = od.fsi_res.near_p1_spec{iC}.subsampled_ppc - od.fsi_res.near_p2_spec{iC}.subsampled_ppc;
                    mean_sts = mean(p1_p2_sts,1);
                    sd_sts = std(p1_p2_sts);
                    mean_ppc = mean(p1_p2_ppc, 1);
                    sd_ppc = std(p1_p2_ppc);
                    hfr_lfr_sts = od.fsi_res.near_hfr_spec{iC}.subsampled_sts - od.fsi_res.near_lfr_spec{iC}.subsampled_sts;
                    hfr_lfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc - od.fsi_res.near_lfr_spec{iC}.subsampled_ppc;
                    
                    % Add "SIG" to filename if any of the hfr_lfr ppc lies outside mean +- 2*sd
                    if sum(hfr_lfr_ppc(x1:end) >= mean_ppc(x1:end) + 2*sd_ppc(x1:end)) ~= 0 & ...
                        sum(hfr_lfr_ppc(x1:end) <= mean_ppc(x1:end) - 2*sd_ppc(x1:end)) ~= 0
                        [~ , f1] = max(hfr_lfr_ppc(x1:end) - (mean_ppc(x1:end) + 2*sd_ppc(x1:end)));
                        [~ , f2] = max((mean_ppc(x1:end) - 2*sd_ppc(x1:end)) - hfr_lfr_ppc(x1:end));
                        foi = od.fsi_res.near_spec{iC}.freqs(x1:end);
                        f1 = foi(f1);
                        f2 = foi(f2);
                        this_prefix = sprintf("_SIGhfr_%.0f_SIGlfr_%.0f",f1,f2);
                        o_prefix = cat(2, o_prefix, convertStringsToChars(this_prefix));
                    elseif sum(hfr_lfr_ppc(x1:end) >= mean_ppc(x1:end) + 2*sd_ppc(x1:end)) ~= 0
                        [~ , f1] = max(hfr_lfr_ppc(x1:end) - (mean_ppc(x1:end) + 2*sd_ppc(x1:end)));
                        foi = od.fsi_res.near_spec{iC}.freqs(x1:end);
                        f1 = foi(f1);
                        this_prefix = sprintf("_SIGhfr_%.0f",f1);
                        o_prefix = cat(2, o_prefix, convertStringsToChars(this_prefix));
                    elseif sum(hfr_lfr_ppc(x1:end) <= mean_ppc(x1:end) - 2*sd_ppc(x1:end)) ~= 0
                        [~ , f2] = max((mean_ppc(x1:end) - 2*sd_ppc(x1:end)) - hfr_lfr_ppc(x1:end));
                        foi = od.fsi_res.near_spec{iC}.freqs(x1:end);
                        f2 = foi(f2);
                        this_prefix = sprintf("_SIGlfr_%.0f",f2);
                        o_prefix = cat(2, o_prefix, convertStringsToChars(this_prefix));
                    end
                              
                    % Plot PPC
                    subplot(3,3,3)
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
                    x1 = find(od.fsi_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    xlabel('Freqs')
                    title('Near Reward PPC');
                    
                    % Plot Control split stuff
                    subplot(3,3,4)
                    plot(od.fsi_res.near_spec{iC}.freqs, mean_sts + 2*sd_sts, '-g')
                    hold on; plot(od.fsi_res.near_spec{iC}.freqs, mean_sts - 2*sd_sts, '-r')
                    plot(od.fsi_res.near_spec{iC}.freqs, hfr_lfr_sts, '--blue')
%                     legend({'P1-P2: Mean + 1*SD', 'P1-P2: Mean - 1*SD', 'HFR - LFR'}, 'FontSize', 8, 'Location', 'southeast')
                    x1 = find(od.fsi_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    ylabel('STS difference')
                    xlabel('Freqs')
                    
                    subplot(3,3,5)
                    plot(od.fsi_res.near_spec{iC}.freqs, mean_ppc + 2*sd_ppc, '-g')
                    hold on; plot(od.fsi_res.near_spec{iC}.freqs, mean_ppc - 2*sd_ppc, '-r')
                    plot(od.fsi_res.near_spec{iC}.freqs, hfr_lfr_ppc, '--blue')
%                     legend({'P1-P2: Mean + 1*SD', 'P1-P2: Mean - 1*SD', 'HFR - LFR'}, 'FontSize', 8, 'Location', 'southeast')
                    x1 = find(od.fsi_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    ylabel('PPC difference')
                    xlabel('Freqs')
                    
                    
                    % Plot histograms
                    subplot(6,6,17)
                    hold on
                    histogram(near_p1_lg_sts_pk - near_p2_lg_sts_pk,lg(1)-lg(2)-0.5:5:lg(2)-lg(1)+0.5,'FaceColor','red','FaceAlpha',0.6);
                    xline(near_hfr_lg_sts_pk - near_lfr_lg_sts_pk, 'blue')
                    ylabel('STS 5-30 Peak dif')
                    subplot(6,6,18)
                    hold on
                    histogram(near_p1_lg_ppc_pk - near_p2_lg_ppc_pk,lg(1)-lg(2)-0.5:5:lg(2)-lg(1)+0.5,'FaceColor','red','FaceAlpha',0.6);
                    xline(near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk, 'blue')
                    ylabel('PPC 5-30 Peak dif')
                    subplot(6,6,23)
                    hold on
                    histogram(near_p1_hg_sts_pk - near_p2_hg_sts_pk,hg(1)-hg(2)-0.5:10:hg(2)-hg(1)+0.5,'FaceColor','green','FaceAlpha',0.6);
                    xline(near_hfr_hg_sts_pk - near_lfr_hg_sts_pk, 'blue')
                    ylabel('STS 30-100 Peak dif')
                    subplot(6,6,24)
                    hold on
                    histogram(near_p1_hg_ppc_pk - near_p2_hg_ppc_pk,hg(1)-hg(2)-0.5:10:hg(2)-hg(1)+0.5,'FaceColor','green','FaceAlpha',0.6);
                    xline(near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk, 'blue')
                    ylabel('PPC 30-100 Peak dif')
                    
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
                    near_p1_text = sprintf("P1 Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_p1_spec_range(1), near_p1_spec_range(2), near_p1_spec_mean);
                    near_p1_text2 = sprintf("MSN P1 spike dist:");
                    near_p1_text3 = sprintf("%d ", od.msn_near_p1_dist(1:end));
                    near_p2_text = sprintf("P2 Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_p2_spec_range(1), near_p2_spec_range(2), near_p2_spec_mean); 
                    near_p2_text2 = sprintf("MSN P2 spike dist: ");
                    near_p2_text3 = sprintf("%d ", od.msn_near_p2_dist(1:end));
                    sim_text = sprintf("LFR-P1 match: %.2f, LFR-P2 match: %.2f", ...
                        length(find(lfr_trials' == p1_trials))/length(p1_trials), ...
                        length(find(lfr_trials' == p2_trials))/length(p2_trials));
                    sim_text2 = sprintf("HFR-P1 match: %.2f, HFR-P2 match: %.2f", ...
                        length(find(hfr_trials' == p1_trials))/length(p1_trials), ...
                        length(find(hfr_trials' == p2_trials))/length(p2_trials));
                    
                    text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.75, near_text2, 'Color', 'blue', 'FontSize', 9);
                    text(0.4, 0.75, near_text3, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.65, near_lfr_text, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.58, near_lfr_text2, 'Color', 'red', 'FontSize', 9);
                    text(0.4, 0.58, near_lfr_text3, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.48, near_hfr_text, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.4, near_hfr_text2, 'Color', 'green', 'FontSize', 9);
                    text(0.4, 0.4, near_hfr_text3, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.3, near_p1_text, 'Color', 'black', 'FontSize', 9);
                    text(0, 0.22, near_p1_text2, 'Color', 'black', 'FontSize', 9);
                    text(0.4, 0.22, near_p1_text3, 'Color', 'black', 'FontSize', 9);
                    text(0, 0.12, near_p2_text, 'Color', 'black', 'FontSize', 9);
                    text(0, 0.04, near_p2_text2, 'Color', 'black', 'FontSize', 9);
                    text(0.4, 0.04, near_p2_text3, 'Color', 'black', 'FontSize', 9);
                    text(0, -0.04, sim_text, 'Color', 'magenta', 'FontSize', 9);
                    text(0, -0.14, sim_text2, 'Color', 'magenta', 'FontSize', 9);
                    axis off;
                    
                    subplot(3,3,8)
                    sts_near_text = sprintf('Correlation between unsampled and subsampled near STS: %.2f', near_subsample_sts_corr);
                    sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
                    sts_lfr_text2 = sprintf('Correlation between unsampled and subsampled LFR STS: %.2f', near_lfr_subsample_sts_corr);
                    sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
                    sts_hfr_text2 = sprintf('Correlation between unsampled and subsampled HFR STS: %.2f', near_hfr_subsample_sts_corr);
                    sts_p1_text = sprintf('Correlation between near and p1 STS: %.2f', near_p1_sts_corr);
                    sts_p1_text2 = sprintf('Correlation between unsampled and subsampled P1 STS: %.2f', near_p1_subsample_sts_corr);
                    sts_p2_text = sprintf('Correlation between near and p2 STS: %.2f', near_p2_sts_corr);
                    sts_p2_text2 = sprintf('Correlation between unsampled and subsampled P2 STS: %.2f', near_p2_subsample_sts_corr);
                    text(0, 0.85, sts_text{1}, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.75, sts_near_text, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.65, sts_text{2}, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.58, sts_lfr_text, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.48, sts_lfr_text2, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.4, sts_text{3}, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.3, sts_hfr_text, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.22, sts_hfr_text2, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.12, sts_p1_text, 'Color', 'magenta', 'FontSize', 9);
                    text(0, 0.04, sts_p1_text2, 'Color', 'magenta', 'FontSize', 9);
                    text(0, -0.04, sts_p2_text, 'Color', 'black', 'FontSize', 9);
                    text(0, -0.14, sts_p2_text2, 'Color', 'black', 'FontSize', 9);
                    axis off;
                    
                    subplot(3,3,9)
                    ppc_near_text = sprintf('Correlation between unsampled and subsampled near PPC: %.2f', near_subsample_ppc_corr);
                    ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
                    ppc_lfr_text2 = sprintf('Correlation between unsampled and subsampled LFR PPC: %.2f', near_lfr_subsample_ppc_corr);
                    ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
                    ppc_hfr_text2 = sprintf('Correlation between unsampled and subsampled HFR PPC: %.2f', near_hfr_subsample_ppc_corr);
                    ppc_p1_text = sprintf('Correlation between near and p1 PPC: %.2f', near_p1_ppc_corr);
                    ppc_p1_text2 = sprintf('Correlation between unsampled and subsampled P1 PPC: %.2f', near_p1_subsample_ppc_corr);
                    ppc_p2_text = sprintf('Correlation between near and p2 PPC: %.2f', near_p2_ppc_corr);
                    ppc_p2_text2 = sprintf('Correlation between unsampled and subsampled P2 PPC: %.2f', near_p2_subsample_ppc_corr);
                    text(0, 0.85, ppc_text{1}, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.75, ppc_near_text, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.65, ppc_text{2}, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.58, ppc_lfr_text, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.48, ppc_lfr_text2, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.4, ppc_text{3}, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.3, ppc_hfr_text, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.22, ppc_hfr_text2, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.12, ppc_p1_text, 'Color', 'magenta', 'FontSize', 9);
                    text(0, 0.04, ppc_p1_text2, 'Color', 'magenta', 'FontSize', 9);
                    text(0, -0.04, ppc_p2_text, 'Color', 'black', 'FontSize', 9);
                    text(0, -0.14, ppc_p2_text2, 'Color', 'black', 'FontSize', 9);
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
            fig = figure('WindowState', 'maximized');
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if isfield(od.msn_res.near_spec{iC},'flag_no_control_split') && ...
                    ~od.msn_res.near_spec{iC}.flag_no_control_split          
                
                % Plot STA
                subplot(3,3,1)
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                flag_leg = true;
                this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.msn_res.near_spec{iC}.spk_count);
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.msn_res.near_lfr_spec{iC}.spk_count);
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.msn_res.near_hfr_spec{iC}.spk_count);
                legend(this_legend, 'Location', 'northwest', 'FontSize', 8)
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
                flag_near_p1_lg_sts_peak = true;
                flag_near_p1_hg_sts_peak = true;
                flag_near_p1_lg_ppc_peak = true;
                flag_near_p1_hg_ppc_peak = true;
                flag_near_p2_lg_sts_peak = true;
                flag_near_p2_hg_sts_peak = true;
                flag_near_p2_lg_ppc_peak = true;
                flag_near_p2_hg_ppc_peak = true;
                near_p1_lg_sts_pk = zeros(1,num_control_splits);
                near_p2_lg_sts_pk = zeros(1,num_control_splits);
                near_p1_hg_sts_pk = zeros(1,num_control_splits);
                near_p2_hg_sts_pk = zeros(1,num_control_splits);
                near_p1_lg_ppc_pk = zeros(1,num_control_splits);
                near_p2_hg_ppc_pk = zeros(1,num_control_splits);
                near_p1_hg_ppc_pk = zeros(1,num_control_splits);
                near_p2_lg_ppc_pk = zeros(1,num_control_splits);
                

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
                
                % Peaks in near_p1_spec
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_p1_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_p1_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p1_spec{iC}.sts(iP,lf:rf));
                    max_val = max(od.msn_res.near_p1_spec{iC}.sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_lg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_lg_sts_peak = false;
                        end
                    else
                        flag_near_p1_lg_sts_peak = false; 
                    end
                end   
                % Find Low Gamma PPC Peak
                lf = find(od.msn_res.near_p1_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_p1_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p1_spec{iC}.ppc(iP,lf:rf));
                    max_val = max(od.msn_res.near_p1_spec{iC}.ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_lg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_lg_ppc_peak = false;
                        end
                    else
                        flag_near_p1_lg_ppc_peak = false; 
                    end
                end
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_p1_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_p1_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p1_spec{iC}.sts(iP,lf:rf));
                    max_val = max(od.msn_res.near_p1_spec{iC}.sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_hg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_hg_sts_peak = false;
                        end
                    else
                        flag_near_p1_hg_sts_peak = false; 
                    end
                end   
                % Find High Gamma PPC Peak
                lf = find(od.msn_res.near_p1_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_p1_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p1_spec{iC}.ppc(iP,lf:rf));
                    max_val = max(od.msn_res.near_p1_spec{iC}.ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p1_hg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p1_hg_ppc_peak = false;
                        end
                    else
                        flag_near_p1_hg_ppc_peak = false; 
                    end
                end
                
                % Peaks in near_p2_spec
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_p2_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_p2_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p2_spec{iC}.sts(iP,lf:rf));
                    max_val = max(od.msn_res.near_p2_spec{iC}.sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_lg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_lg_sts_peak = false;
                        end
                    else
                        flag_near_p2_lg_sts_peak = false; 
                    end
                end   
                % Find Low Gamma PPC Peak
                lf = find(od.msn_res.near_p2_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_p2_spec{iC}.freqs <= lg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p2_spec{iC}.ppc(iP,lf:rf));
                    max_val = max(od.msn_res.near_p2_spec{iC}.ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_lg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_lg_ppc_peak = false;
                        end
                    else
                        flag_near_p2_lg_ppc_peak = false; 
                    end
                end
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_p2_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_p2_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p2_spec{iC}.sts(iP,lf:rf));
                    max_val = max(od.msn_res.near_p2_spec{iC}.sts(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.sts(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_hg_sts_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_hg_sts_peak = false;
                        end
                    else
                        flag_near_p2_hg_sts_peak = false; 
                    end
                end   
                % Find High Gamma PPC Peak
                lf = find(od.msn_res.near_p2_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_p2_spec{iC}.freqs <= hg(2), 1, 'last');
                for iP = 1:num_control_splits
                    pks = findpeaks(od.msn_res.near_p2_spec{iC}.ppc(iP,lf:rf));
                    max_val = max(od.msn_res.near_p2_spec{iC}.ppc(iP,:));
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.ppc(iP,lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            near_p2_hg_ppc_pk(iP) = pks.loc(loc1);
                        else
                           flag_near_p2_hg_ppc_peak = false;
                        end
                    else
                        flag_near_p2_hg_ppc_peak = false; 
                    end
                end
                
                all_sts_peaks_in_lg = flag_near_lg_sts_peak && ...
                    flag_near_hfr_lg_sts_peak && ...
                    flag_near_lfr_lg_sts_peak && ...
                    flag_near_p2_lg_sts_peak && ...
                    flag_near_p1_lg_sts_peak;
                all_sts_peaks_in_hg = flag_near_hg_sts_peak && ...
                    flag_near_hfr_hg_sts_peak && ...
                    flag_near_lfr_hg_sts_peak && ...
                    flag_near_p2_hg_sts_peak && ...
                    flag_near_p1_hg_sts_peak;
                all_ppc_peaks_in_lg = flag_near_lg_ppc_peak && ...
                    flag_near_hfr_lg_ppc_peak && ...
                    flag_near_lfr_lg_ppc_peak && ...
                    flag_near_p2_lg_ppc_peak && ...
                    flag_near_p1_lg_ppc_peak;                              
                all_ppc_peaks_in_hg = flag_near_hg_ppc_peak && ...
                    flag_near_hfr_hg_ppc_peak && ...
                    flag_near_lfr_hg_ppc_peak && ...
                    flag_near_p2_hg_ppc_peak && ...
                    flag_near_p1_hg_ppc_peak;

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
                near_lfr_sts_corr = corrcoef(od.msn_res.near_spec{iC}.sts_vals,od.msn_res.near_lfr_spec{iC}.sts_vals);
                near_lfr_sts_corr = near_lfr_sts_corr(1,2);
                near_hfr_sts_corr = corrcoef(od.msn_res.near_spec{iC}.sts_vals,od.msn_res.near_hfr_spec{iC}.sts_vals);
                near_hfr_sts_corr = near_hfr_sts_corr(1,2);
                near_lfr_ppc_corr = corrcoef(od.msn_res.near_spec{iC}.ppc,od.msn_res.near_lfr_spec{iC}.ppc);
                near_lfr_ppc_corr = near_lfr_ppc_corr(1,2);
                near_hfr_ppc_corr = corrcoef(od.msn_res.near_spec{iC}.ppc,od.msn_res.near_hfr_spec{iC}.ppc);
                near_hfr_ppc_corr = near_hfr_ppc_corr(1,2);
                 
                p1_trials = od.msn_res.near_spec{iC}.valid_splits;    
                p2_trials = ~p1_trials;
                p1_trials = p1_trials & repmat(nz_trials',100,1);
                p2_trials = p2_trials & repmat(nz_trials',100,1);
                temp_p1_spk = zeros(100,1);
                temp_p1_dur = zeros(100,1);
                temp_p2_spk = zeros(100,1);
                temp_p2_dur = zeros(100,1);
                for iP = 1:100
                    temp_p1_spk(iP) = sum(od.msn_res.near_spec{iC}.trial_spk_count(p1_trials(iP,:)));
                    temp_p1_dur(iP) = sum(od.msn_res.near_spec{iC}.trial_spk_count(p1_trials(iP,:))./od.msn_res.near_spec{iC}.mfr(p1_trials(iP,:))');
                    temp_p2_spk(iP) = sum(od.msn_res.near_spec{iC}.trial_spk_count(p2_trials(iP,:)));
                    temp_p2_dur(iP) = sum(od.msn_res.near_spec{iC}.trial_spk_count(p2_trials(iP,:))./od.msn_res.near_spec{iC}.mfr(p2_trials(iP,:))');
                end
                temp_p1_mfr = temp_p1_spk./temp_p1_dur;
                temp_p2_mfr = temp_p2_spk./temp_p2_dur;
                near_p1_spec_mean = mean(temp_p1_mfr);
                near_p1_spec_range = [min(temp_p1_mfr) max(temp_p1_mfr)];
                near_p2_spec_mean = mean(temp_p2_mfr);
                near_p2_spec_range = [min(temp_p2_mfr) max(temp_p2_mfr)];
                near_p1_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,mean(od.msn_res.near_p1_spec{iC}.sts,1));
                near_p1_sts_corr = near_p1_sts_corr(1,2);
                near_p2_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,mean(od.msn_res.near_p2_spec{iC}.sts,1));
                near_p2_sts_corr = near_p2_sts_corr(1,2);
                near_p1_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,mean(od.msn_res.near_p1_spec{iC}.ppc,1));
                near_p1_ppc_corr = near_p1_ppc_corr(1,2);
                near_p2_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,mean(od.msn_res.near_p2_spec{iC}.ppc,1));
                near_p2_ppc_corr = near_p2_ppc_corr(1,2);
              
                % if peaks exist in both low gamma and high gamma range
                if all_sts_peaks_in_lg && all_ppc_peaks_in_lg && ...
                    all_sts_peaks_in_hg && all_ppc_peaks_in_hg ...
                    
                    % Plot STS
                    subplot(3,3,2)
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
                    x1 = find(od.msn_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    xlabel('Freqs')
                    title('Near Reward STS');
                    
                    p1_p2_sts = od.msn_res.near_p1_spec{iC}.sts - od.msn_res.near_p2_spec{iC}.sts;
                    p1_p2_ppc = od.msn_res.near_p1_spec{iC}.ppc - od.msn_res.near_p2_spec{iC}.ppc;
                    mean_sts = mean(p1_p2_sts,1);
                    sd_sts = std(p1_p2_sts);
                    mean_ppc = mean(p1_p2_ppc, 1);
                    sd_ppc = std(p1_p2_ppc);
                    hfr_lfr_sts = od.msn_res.near_hfr_spec{iC}.sts_vals - od.msn_res.near_lfr_spec{iC}.sts_vals;
                    hfr_lfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc - od.msn_res.near_lfr_spec{iC}.ppc;
                    
                    % Add "SIG" to filename if any of the hfr_lfr ppc lies outside mean +- 2*sd
                    if sum(hfr_lfr_ppc(x1:end)' >= mean_ppc(x1:end) + 2*sd_ppc(x1:end)) ~= 0 & ...
                        sum(hfr_lfr_ppc(x1:end)' <= mean_ppc(x1:end) - 2*sd_ppc(x1:end)) ~= 0
                        [~ , f1] = max(hfr_lfr_ppc(x1:end)' - (mean_ppc(x1:end) + 2*sd_ppc(x1:end)));
                        [~ , f2] = max((mean_ppc(x1:end) - 2*sd_ppc(x1:end)) - hfr_lfr_ppc(x1:end)');
                        foi = od.msn_res.near_spec{iC}.freqs(x1:end);
                        f1 = foi(f1);
                        f2 = foi(f2);
                        this_prefix = sprintf("_SIGhfr_%.0f_SIGlfr_%.0f",f1,f2);
                        o_prefix = cat(2, o_prefix, convertStringsToChars(this_prefix));
                    elseif sum(hfr_lfr_ppc(x1:end)' >= mean_ppc(x1:end) + 2*sd_ppc(x1:end)) ~= 0
                        [~ , f1] = max(hfr_lfr_ppc(x1:end)' - (mean_ppc(x1:end) + 2*sd_ppc(x1:end)));
                        foi = od.msn_res.near_spec{iC}.freqs(x1:end);
                        f1 = foi(f1);
                        this_prefix = sprintf("_SIGhfr_%.0f",f1);
                        o_prefix = cat(2, o_prefix, convertStringsToChars(this_prefix));
                    elseif sum(hfr_lfr_ppc(x1:end)' <= mean_ppc(x1:end) - 2*sd_ppc(x1:end)) ~= 0
                        [~ , f2] = max((mean_ppc(x1:end) - 2*sd_ppc(x1:end)) - hfr_lfr_ppc(x1:end)');
                        foi = od.msn_res.near_spec{iC}.freqs(x1:end);
                        f2 = foi(f2);
                        this_prefix = sprintf("_SIGlfr_%.0f",f2);
                        o_prefix = cat(2, o_prefix, convertStringsToChars(this_prefix));
                    end
                                   
                    % Plot PPC
                    subplot(3,3,3)
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
                    x1 = find(od.msn_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    xlabel('Freqs')
                    title('Near Reward PPC');
                    
                    % Plot Control split stuff
                    subplot(3,3,4)
                    plot(od.msn_res.near_spec{iC}.freqs, mean_sts + 2*sd_sts, '-g')
                    hold on; plot(od.msn_res.near_spec{iC}.freqs, mean_sts - 2*sd_sts, '-r')
                    plot(od.msn_res.near_spec{iC}.freqs, hfr_lfr_sts, '--blue')
%                     legend({'P1-P2: Mean + 1*SD', 'P1-P2: Mean - 1*SD', 'HFR - LFR'}, 'FontSize', 8, 'Location', 'southeast')
                    x1 = find(od.msn_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    ylabel('STS difference')
                    xlabel('Freqs')
                    
                    subplot(3,3,5)
                    plot(od.msn_res.near_spec{iC}.freqs, mean_ppc + 2*sd_ppc, '-g')
                    hold on; plot(od.msn_res.near_spec{iC}.freqs, mean_ppc - 2*sd_ppc, '-r')
                    plot(od.msn_res.near_spec{iC}.freqs, hfr_lfr_ppc, '--blue')
%                     legend({'P1-P2: Mean + 1*SD', 'P1-P2: Mean - 1*SD', 'HFR - LFR'}, 'FontSize', 8, 'Location', 'southeast')
                    x1 = find(od.msn_res.near_spec{iC}.freqs >= 5, 1, 'first');
                    x2 = 100;
                    xlim([x1 x2])
                    ylabel('PPC difference')
                    xlabel('Freqs')
                    
                    % Plot histograms
                    subplot(6,6,17)
                    hold on
                    histogram(near_p1_lg_sts_pk - near_p2_lg_sts_pk,lg(1)-lg(2)-0.5:5:lg(2)-lg(1)+0.5,'FaceColor','red','FaceAlpha',0.6);
                    xline(near_hfr_lg_sts_pk - near_lfr_lg_sts_pk, 'blue')
                    ylabel('STS 5-30 Peak dif')
                    subplot(6,6,18)
                    hold on
                    histogram(near_p1_lg_ppc_pk - near_p2_lg_ppc_pk,lg(1)-lg(2)-0.5:5:lg(2)-lg(1)+0.5,'FaceColor','red','FaceAlpha',0.6);
                    xline(near_hfr_lg_ppc_pk - near_lfr_lg_ppc_pk, 'blue')
                    ylabel('PPC 5-30 Peak dif')
                    subplot(6,6,23)
                    hold on
                    histogram(near_p1_hg_sts_pk - near_p2_hg_sts_pk,hg(1)-hg(2)-0.5:10:hg(2)-hg(1)+0.5,'FaceColor','green','FaceAlpha',0.6);
                    xline(near_hfr_hg_sts_pk - near_lfr_hg_sts_pk, 'blue')
                    ylabel('STS 30-100 Peak dif')
                    subplot(6,6,24)
                    hold on
                    histogram(near_p1_hg_ppc_pk - near_p2_hg_ppc_pk,hg(1)-hg(2)-0.5:10:hg(2)-hg(1)+0.5,'FaceColor','green','FaceAlpha',0.6);
                    xline(near_hfr_hg_ppc_pk - near_lfr_hg_ppc_pk, 'blue')
                    ylabel('PPC 30-100 Peak dif')
                    
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
                    near_p1_text = sprintf("P1 Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_p1_spec_range(1), near_p1_spec_range(2), near_p1_spec_mean);
                    near_p1_text2 = sprintf("MSN P1 spike dist:");
                    near_p1_text3 = sprintf("%d ", od.msn_near_p1_dist(1:end));
                    near_p2_text = sprintf("P2 Trials firing rate range %.2f Hz - %.2f Hz, mean fr: %.2f Hz", ...
                        near_p2_spec_range(1), near_p2_spec_range(2), near_p2_spec_mean); 
                    near_p2_text2 = sprintf("MSN P2 spike dist: ");
                    near_p2_text3 = sprintf("%d ", od.msn_near_p2_dist(1:end));
                    sim_text = sprintf("LFR-P1 match: %.2f, LFR-P2 match: %.2f", ...
                        length(find(lfr_trials' == p1_trials))/length(p1_trials), ...
                        length(find(lfr_trials' == p2_trials))/length(p2_trials));
                    sim_text2 = sprintf("HFR-P1 match: %.2f, HFR-P2 match: %.2f", ...
                        length(find(hfr_trials' == p1_trials))/length(p1_trials), ...
                        length(find(hfr_trials' == p2_trials))/length(p2_trials));
                    
                    text(0, 0.8, near_text, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.75, near_text2, 'Color', 'blue', 'FontSize', 9);
                    text(0.4, 0.75, near_text3, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.65, near_lfr_text, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.58, near_lfr_text2, 'Color', 'red', 'FontSize', 9);
                    text(0.4, 0.58, near_lfr_text3, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.48, near_hfr_text, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.4, near_hfr_text2, 'Color', 'green', 'FontSize', 9);
                    text(0.4, 0.4, near_hfr_text3, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.3, near_p1_text, 'Color', 'black', 'FontSize', 9);
                    text(0, 0.22, near_p1_text2, 'Color', 'black', 'FontSize', 9);
                    text(0.4, 0.22, near_p1_text3, 'Color', 'black', 'FontSize', 9);
                    text(0, 0.12, near_p2_text, 'Color', 'black', 'FontSize', 9);
                    text(0, 0.04, near_p2_text2, 'Color', 'black', 'FontSize', 9);
                    text(0.4, 0.04, near_p2_text3, 'Color', 'black', 'FontSize', 9);
                    text(0, -0.04, sim_text, 'Color', 'magenta', 'FontSize', 9);
                    text(0, -0.14, sim_text2, 'Color', 'magenta', 'FontSize', 9);
                    axis off;
                    
                    subplot(3,3,8)
                    sts_lfr_text = sprintf('Correlation between near and lfr STS: %.2f', near_lfr_sts_corr);
                    sts_hfr_text = sprintf('Correlation between near and hfr STS: %.2f', near_hfr_sts_corr);
                    sts_p1_text = sprintf('Correlation between near and p1 STS: %.2f', near_p1_sts_corr);
                    sts_p2_text = sprintf('Correlation between near and p2 STS: %.2f', near_p2_sts_corr);
                    text(0, 1, sts_text{1}, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.9, sts_text{2}, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.8, sts_lfr_text, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.7, sts_text{3}, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.6, sts_hfr_text, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.4, sts_p1_text, 'Color', 'magenta', 'FontSize', 9);
                    text(0, 0.2, sts_p2_text, 'Color', 'black', 'FontSize', 9);
                    axis off;
                    
                    subplot(3,3,9)
                    ppc_lfr_text = sprintf('Correlation between near and lfr PPC: %.2f', near_lfr_ppc_corr);
                    ppc_hfr_text = sprintf('Correlation between near and hfr PPC: %.2f', near_hfr_ppc_corr);
                    ppc_p1_text = sprintf('Correlation between near and p1 PPC: %.2f', near_p1_ppc_corr);
                    ppc_p2_text = sprintf('Correlation between near and p2 PPC: %.2f', near_p2_ppc_corr);
                    text(0, 1, ppc_text{1}, 'Color', 'blue', 'FontSize', 9);
                    text(0, 0.9, ppc_text{2}, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.8, ppc_lfr_text, 'Color', 'red', 'FontSize', 9);
                    text(0, 0.7, ppc_text{3}, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.6, ppc_hfr_text, 'Color', 'green', 'FontSize', 9);
                    text(0, 0.4, ppc_p1_text, 'Color', 'magenta', 'FontSize', 9);
                    text(0, 0.2, ppc_p2_text, 'Color', 'black', 'FontSize', 9);
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