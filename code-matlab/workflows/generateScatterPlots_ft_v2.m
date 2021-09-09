fsi_near_sts_lg_peaks = [];
fsi_near_sts_hg_peaks = [];
fsi_near_lfr_hg_sts_peaks = [];
fsi_near_lfr_lg_sts_peaks = [];
fsi_near_hfr_hg_sts_peaks = [];
fsi_near_hfr_lg_sts_peaks = [];
fsi_near_ppc_lg_peaks = [];
fsi_near_ppc_hg_peaks = [];
fsi_near_lfr_hg_ppc_peaks = [];
fsi_near_lfr_lg_ppc_peaks = [];
fsi_near_hfr_hg_ppc_peaks = [];
fsi_near_hfr_lg_ppc_peaks = [];
fsi_near_lfr_sts_corrs = [];
fsi_near_hfr_sts_corrs = [];
fsi_near_lfr_ppc_corrs = [];
fsi_near_hfr_ppc_corrs = [];
fsi_near_frs = [];
fsi_near_lfr_frs = [];
fsi_near_hfr_frs = [];
fsi_near_p1_hg_sts_peaks = [];
fsi_near_p1_lg_sts_peaks = [];
fsi_near_p2_hg_sts_peaks = [];
fsi_near_p2_lg_sts_peaks = [];
fsi_near_p1_hg_ppc_peaks = [];
fsi_near_p1_lg_ppc_peaks = [];
fsi_near_p2_hg_ppc_peaks = [];
fsi_near_p2_lg_ppc_peaks = [];
fsi_near_p1_sts_corrs = [];
fsi_near_p2_sts_corrs = [];
fsi_near_p1_ppc_corrs = [];
fsi_near_p2_ppc_corrs = [];
fsi_near_p1_frs = [];
fsi_near_p2_frs = [];

msn_near_sts_lg_peaks = [];
msn_near_sts_hg_peaks = [];
msn_near_lfr_hg_sts_peaks = [];
msn_near_lfr_lg_sts_peaks = [];
msn_near_hfr_hg_sts_peaks = [];
msn_near_hfr_lg_sts_peaks = [];
msn_near_ppc_lg_peaks = [];
msn_near_ppc_hg_peaks = [];
msn_near_lfr_hg_ppc_peaks = [];
msn_near_lfr_lg_ppc_peaks = [];
msn_near_hfr_hg_ppc_peaks = [];
msn_near_hfr_lg_ppc_peaks = [];
msn_near_lfr_sts_corrs = [];
msn_near_hfr_sts_corrs = [];
msn_near_lfr_ppc_corrs = [];
msn_near_hfr_ppc_corrs = [];
msn_near_frs = [];
msn_near_lfr_frs = [];
msn_near_hfr_frs = [];
msn_near_p1_hg_sts_peaks = [];
msn_near_p1_lg_sts_peaks = [];
msn_near_p2_hg_sts_peaks = [];
msn_near_p2_lg_sts_peaks = [];
msn_near_p1_hg_ppc_peaks = [];
msn_near_p1_lg_ppc_peaks = [];
msn_near_p2_hg_ppc_peaks = [];
msn_near_p2_lg_ppc_peaks = [];
msn_near_p1_sts_corrs = [];
msn_near_p2_sts_corrs = [];
msn_near_p1_ppc_corrs = [];
msn_near_p2_ppc_corrs = [];
msn_near_p1_frs = [];
msn_near_p2_frs = [];

sessions_with_msn_ppc_peak_ge_fsi = [];

cd('D:\RandomVstrAnalysis\ft_results');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
lg = [30,80];
hg = [30,80];
pk_thresh = -1;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_near_sts_lg_peak_vals = [];
        fsi_near_sts_hg_peak_vals = [];
        fsi_near_lfr_hg_sts_peak_vals = [];
        fsi_near_lfr_lg_sts_peak_vals = [];
        fsi_near_hfr_hg_sts_peak_vals = [];
        fsi_near_hfr_lg_sts_peak_vals = [];
        fsi_near_ppc_lg_peak_vals = [];
        fsi_near_ppc_hg_peak_vals = [];
        fsi_near_lfr_hg_ppc_peak_vals = [];
        fsi_near_lfr_lg_ppc_peak_vals = [];
        fsi_near_hfr_hg_ppc_peak_vals = [];
        fsi_near_hfr_lg_ppc_peak_vals = [];
        fsi_near_p1_hg_sts_peak_vals = [];
        fsi_near_p1_lg_sts_peak_vals = [];
        fsi_near_p2_hg_sts_peak_vals = [];
        fsi_near_p2_lg_sts_peak_vals = [];
        fsi_near_p1_hg_ppc_peak_vals = [];
        fsi_near_p1_lg_ppc_peak_vals = [];
        fsi_near_p2_hg_ppc_peak_vals = [];
        fsi_near_p2_lg_ppc_peak_vals = [];
        for iC = 1:length(fsi_labels)
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
                near_p1_subsample_ppc_corr = corrcoef(od.fsi_res.near_p1_spec{iC}.ppc ,od.fsi_res.near_p1_spec{iC}.subsampled_ppc);
                near_p1_subsample_ppc_corr = near_p1_subsample_ppc_corr(1,2);
                near_p2_subsample_ppc_corr = corrcoef(od.fsi_res.near_p2_spec{iC}.ppc ,od.fsi_res.near_p2_spec{iC}.subsampled_ppc);
                near_p2_subsample_ppc_corr = near_p2_subsample_ppc_corr(1,2);
                near_subsample_sts_corr = corrcoef(od.fsi_res.near_spec{iC}.sts_vals ,od.fsi_res.near_spec{iC}.subsampled_sts);
                near_subsample_sts_corr = near_subsample_sts_corr(1,2);
                near_lfr_subsample_sts_corr = corrcoef(od.fsi_res.near_lfr_spec{iC}.sts_vals ,od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                near_lfr_subsample_sts_corr = near_lfr_subsample_sts_corr(1,2);
                near_hfr_subsample_sts_corr = corrcoef(od.fsi_res.near_hfr_spec{iC}.sts_vals ,od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                near_hfr_subsample_sts_corr = near_hfr_subsample_sts_corr(1,2);
                near_p1_subsample_sts_corr = corrcoef(od.fsi_res.near_p1_spec{iC}.sts_vals ,od.fsi_res.near_p1_spec{iC}.subsampled_sts);
                near_p1_subsample_sts_corr = near_p1_subsample_sts_corr(1,2);
                near_p2_subsample_sts_corr = corrcoef(od.fsi_res.near_p2_spec{iC}.sts_vals ,od.fsi_res.near_p2_spec{iC}.subsampled_sts);
                near_p2_subsample_sts_corr = near_p2_subsample_sts_corr(1,2);
                
                % Skip if any subsampled correlations are less than 0.9
                if near_subsample_ppc_corr < 0.9 || near_lfr_subsample_ppc_corr < 0.9 || near_hfr_subsample_ppc_corr < 0.9 || ...
                    near_subsample_sts_corr < 0.9 || near_lfr_subsample_sts_corr < 0.9 || near_hfr_subsample_sts_corr < 0.9 || ...
                    near_p1_subsample_ppc_corr < 0.9 || near_p2_subsample_ppc_corr < 0.9 || ...
                    near_p1_subsample_sts_corr < 0.9 || near_p2_subsample_sts_corr < 0.9
                    close all;
                    continue
                end
                
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
                pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p1_lg_sts_peak = false;
                    end
                else
                    flag_near_p1_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p1_lg_ppc_peak = false;
                    end
                else
                    flag_near_p1_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_p1_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_p1_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p1_hg_sts_peak = false;
                    end
                else
                    flag_near_p1_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p1_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p1_hg_ppc_peak = false;
                    end
                else
                    flag_near_p1_hg_ppc_peak = false; 
                end
                
                % Peaks in near_p2_spec
                % Find Low Gamma STS Peak
                lf = find(od.fsi_res.near_p2_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_p2_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p2_lg_sts_peak = false;
                    end
                else
                    flag_near_p2_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p2_lg_ppc_peak = false;
                    end
                else
                    flag_near_p2_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.fsi_res.near_p2_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_p2_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_sts(lf:rf));
                max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_sts(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p2_hg_sts_peak = false;
                    end
                else
                    flag_near_p2_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(lf:rf));
                max_val = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.fsi_res.near_p2_spec{iC}.subsampled_ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p2_hg_ppc_peak = false;
                    end
                else
                    flag_near_p2_hg_ppc_peak = false; 
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
                near_lfr_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,od.fsi_res.near_lfr_spec{iC}.subsampled_sts);
                near_lfr_sts_corr = near_lfr_sts_corr(1,2);
                near_hfr_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,od.fsi_res.near_hfr_spec{iC}.subsampled_sts);
                near_hfr_sts_corr = near_hfr_sts_corr(1,2);
                near_lfr_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,od.fsi_res.near_lfr_spec{iC}.subsampled_ppc);
                near_lfr_ppc_corr = near_lfr_ppc_corr(1,2);
                near_hfr_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,od.fsi_res.near_hfr_spec{iC}.subsampled_ppc);
                near_hfr_ppc_corr = near_hfr_ppc_corr(1,2);
                 
                p1_trials = od.fsi_res.near_spec{iC}.p1_trials;     
                p2_trials = od.fsi_res.near_spec{iC}.p2_trials;  
                p1_trials = p1_trials & nz_trials';
                p2_trials = p2_trials & nz_trials';
                near_p1_spec_mean = sum(od.fsi_res.near_spec{iC}.trial_spk_count(p1_trials))/sum(od.fsi_res.near_spec{iC}.trial_spk_count(p1_trials)'./od.fsi_res.near_spec{iC}.mfr(p1_trials));
                near_p1_spec_range = [min(od.fsi_res.near_spec{iC}.mfr(p1_trials)), max(od.fsi_res.near_spec{iC}.mfr(p1_trials))];
                near_p2_spec_mean = sum(od.fsi_res.near_spec{iC}.trial_spk_count(p2_trials))/sum(od.fsi_res.near_spec{iC}.trial_spk_count(p2_trials)'./od.fsi_res.near_spec{iC}.mfr(p2_trials));
                near_p2_spec_range = [min(od.fsi_res.near_spec{iC}.mfr(p2_trials)), max(od.fsi_res.near_spec{iC}.mfr(p2_trials))];
                near_p1_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,od.fsi_res.near_p1_spec{iC}.subsampled_sts);
                near_p1_sts_corr = near_p1_sts_corr(1,2);
                near_p2_sts_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_sts ,od.fsi_res.near_p2_spec{iC}.subsampled_sts);
                near_p2_sts_corr = near_p2_sts_corr(1,2);
                near_p1_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,od.fsi_res.near_p1_spec{iC}.subsampled_ppc);
                near_p1_ppc_corr = near_p1_ppc_corr(1,2);
                near_p2_ppc_corr  = corrcoef(od.fsi_res.near_spec{iC}.subsampled_ppc ,od.fsi_res.near_p2_spec{iC}.subsampled_ppc);
                near_p2_ppc_corr = near_p2_ppc_corr(1,2);
              
                % if peaks exist in both low gamma and high gamma range
                if all_sts_peaks_in_lg && all_ppc_peaks_in_lg && ...
                    all_sts_peaks_in_hg && all_ppc_peaks_in_hg ...
                    
                    near_lg_sts_pk_val = od.fsi_res.near_spec{iC}.subsampled_sts(lg(1)+near_lg_sts_pk-1);
                    near_hg_sts_pk_val = od.fsi_res.near_spec{iC}.subsampled_sts(hg(1)+near_hg_sts_pk-1);
                    near_lfr_lg_sts_pk_val = od.fsi_res.near_lfr_spec{iC}.subsampled_sts(lg(1)+near_lfr_lg_sts_pk-1);
                    near_lfr_hg_sts_pk_val = od.fsi_res.near_lfr_spec{iC}.subsampled_sts(hg(1)+near_lfr_hg_sts_pk-1);
                    near_hfr_lg_sts_pk_val = od.fsi_res.near_hfr_spec{iC}.subsampled_sts(lg(1)+near_hfr_lg_sts_pk-1);
                    near_hfr_hg_sts_pk_val = od.fsi_res.near_hfr_spec{iC}.subsampled_sts(hg(1)+near_hfr_hg_sts_pk-1);
                    near_lg_ppc_pk_val = od.fsi_res.near_spec{iC}.subsampled_ppc(lg(1)+near_lg_ppc_pk-1);
                    near_hg_ppc_pk_val = od.fsi_res.near_spec{iC}.subsampled_ppc(hg(1)+near_hg_ppc_pk-1);
                    near_lfr_lg_ppc_pk_val = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(lg(1)+near_lfr_lg_ppc_pk-1);
                    near_lfr_hg_ppc_pk_val = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(hg(1)+near_lfr_hg_ppc_pk-1);
                    near_hfr_lg_ppc_pk_val = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(lg(1)+near_hfr_lg_ppc_pk-1);
                    near_hfr_hg_ppc_pk_val = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(hg(1)+near_hfr_hg_ppc_pk-1);
                    near_p1_lg_sts_pk_val = od.fsi_res.near_p1_spec{iC}.subsampled_sts(lg(1)+near_p1_lg_sts_pk-1);
                    near_p1_hg_sts_pk_val = od.fsi_res.near_p1_spec{iC}.subsampled_sts(hg(1)+near_p1_hg_sts_pk-1);
                    near_p2_lg_sts_pk_val = od.fsi_res.near_p2_spec{iC}.subsampled_sts(lg(1)+near_p2_lg_sts_pk-1);
                    near_p2_hg_sts_pk_val = od.fsi_res.near_p2_spec{iC}.subsampled_sts(hg(1)+near_p2_hg_sts_pk-1);
                    near_p1_lg_ppc_pk_val = od.fsi_res.near_p1_spec{iC}.subsampled_ppc(lg(1)+near_p1_lg_ppc_pk-1);
                    near_p1_hg_ppc_pk_val = od.fsi_res.near_p1_spec{iC}.subsampled_ppc(hg(1)+near_p1_hg_ppc_pk-1);
                    near_p2_lg_ppc_pk_val = od.fsi_res.near_p2_spec{iC}.subsampled_ppc(lg(1)+near_p2_lg_ppc_pk-1);
                    near_p2_hg_ppc_pk_val = od.fsi_res.near_p2_spec{iC}.subsampled_ppc(hg(1)+near_p2_hg_ppc_pk-1);
                    
                    near_lg_sts_pk = od.fsi_res.near_spec{iC}.freqs(lg(1)+near_lg_sts_pk-1);
                    near_hg_sts_pk = od.fsi_res.near_spec{iC}.freqs(hg(1)+near_hg_sts_pk-1);
                    near_lfr_lg_sts_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_sts_pk-1);
                    near_lfr_hg_sts_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_sts_pk-1);
                    near_hfr_lg_sts_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_sts_pk-1);
                    near_hfr_hg_sts_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_sts_pk-1);
                    near_lg_ppc_pk = od.fsi_res.near_spec{iC}.freqs(lg(1)+near_lg_ppc_pk-1);
                    near_hg_ppc_pk = od.fsi_res.near_spec{iC}.freqs(hg(1)+near_hg_ppc_pk-1);
                    near_lfr_lg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
                    near_lfr_hg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
                    near_hfr_lg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
                    near_hfr_hg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);
                    near_p1_lg_sts_pk = od.fsi_res.near_p1_spec{iC}.freqs(lg(1)+near_p1_lg_sts_pk-1);
                    near_p1_hg_sts_pk = od.fsi_res.near_p1_spec{iC}.freqs(hg(1)+near_p1_hg_sts_pk-1);
                    near_p2_lg_sts_pk = od.fsi_res.near_p2_spec{iC}.freqs(lg(1)+near_p2_lg_sts_pk-1);
                    near_p2_hg_sts_pk = od.fsi_res.near_p2_spec{iC}.freqs(hg(1)+near_p2_hg_sts_pk-1);
                    near_p1_lg_ppc_pk = od.fsi_res.near_p1_spec{iC}.freqs(lg(1)+near_p1_lg_ppc_pk-1);
                    near_p1_hg_ppc_pk = od.fsi_res.near_p1_spec{iC}.freqs(hg(1)+near_p1_hg_ppc_pk-1);
                    near_p2_lg_ppc_pk = od.fsi_res.near_p2_spec{iC}.freqs(lg(1)+near_p2_lg_ppc_pk-1);
                    near_p2_hg_ppc_pk = od.fsi_res.near_p2_spec{iC}.freqs(hg(1)+near_p2_hg_ppc_pk-1);
                    
                    fsi_near_sts_lg_peaks = [fsi_near_sts_lg_peaks, near_lg_sts_pk];
                    fsi_near_sts_hg_peaks = [fsi_near_sts_hg_peaks, near_hg_sts_pk];
                    fsi_near_lfr_hg_sts_peaks = [fsi_near_lfr_hg_sts_peaks, near_lfr_hg_sts_pk];
                    fsi_near_lfr_lg_sts_peaks = [fsi_near_lfr_lg_sts_peaks, near_lfr_lg_sts_pk];
                    fsi_near_hfr_hg_sts_peaks = [fsi_near_hfr_hg_sts_peaks, near_hfr_hg_sts_pk];
                    fsi_near_hfr_lg_sts_peaks = [fsi_near_hfr_lg_sts_peaks, near_hfr_lg_sts_pk];
                    fsi_near_ppc_lg_peaks = [fsi_near_ppc_lg_peaks, near_lg_ppc_pk];
                    fsi_near_ppc_hg_peaks = [fsi_near_ppc_hg_peaks, near_hg_ppc_pk];
                    fsi_near_lfr_hg_ppc_peaks = [fsi_near_lfr_hg_ppc_peaks, near_lfr_hg_ppc_pk];
                    fsi_near_lfr_lg_ppc_peaks = [fsi_near_lfr_lg_ppc_peaks, near_lfr_lg_ppc_pk];
                    fsi_near_hfr_hg_ppc_peaks = [fsi_near_hfr_hg_ppc_peaks, near_hfr_hg_ppc_pk];
                    fsi_near_hfr_lg_ppc_peaks = [fsi_near_hfr_lg_ppc_peaks, near_hfr_lg_ppc_pk];
                    fsi_near_p1_hg_sts_peaks = [fsi_near_p1_hg_sts_peaks, near_p1_hg_sts_pk];
                    fsi_near_p1_lg_sts_peaks = [fsi_near_p1_lg_sts_peaks, near_p1_lg_sts_pk];
                    fsi_near_p2_hg_sts_peaks = [fsi_near_p2_hg_sts_peaks, near_p2_hg_sts_pk];
                    fsi_near_p2_lg_sts_peaks = [fsi_near_p2_lg_sts_peaks, near_p2_lg_sts_pk];
                    fsi_near_p1_hg_ppc_peaks = [fsi_near_p1_hg_ppc_peaks, near_p1_hg_ppc_pk];
                    fsi_near_p1_lg_ppc_peaks = [fsi_near_p1_lg_ppc_peaks, near_p1_lg_ppc_pk];
                    fsi_near_p2_hg_ppc_peaks = [fsi_near_p2_hg_ppc_peaks, near_p2_hg_ppc_pk];
                    fsi_near_p2_lg_ppc_peaks = [fsi_near_p2_lg_ppc_peaks, near_p2_lg_ppc_pk];
                    
                    fsi_near_sts_lg_peak_vals = [fsi_near_sts_lg_peak_vals, near_lg_sts_pk_val];
                    fsi_near_sts_hg_peak_vals = [fsi_near_sts_hg_peak_vals, near_hg_sts_pk_val];
                    fsi_near_lfr_hg_sts_peak_vals = [fsi_near_lfr_hg_sts_peak_vals, near_lfr_hg_sts_pk_val];
                    fsi_near_lfr_lg_sts_peak_vals = [fsi_near_lfr_lg_sts_peak_vals, near_lfr_lg_sts_pk_val];
                    fsi_near_hfr_hg_sts_peak_vals = [fsi_near_hfr_hg_sts_peak_vals, near_hfr_hg_sts_pk_val];
                    fsi_near_hfr_lg_sts_peak_vals = [fsi_near_hfr_lg_sts_peak_vals, near_hfr_lg_sts_pk_val];
                    fsi_near_ppc_lg_peak_vals = [fsi_near_ppc_lg_peak_vals, near_lg_ppc_pk_val];
                    fsi_near_ppc_hg_peak_vals = [fsi_near_ppc_hg_peak_vals, near_hg_ppc_pk_val];
                    fsi_near_lfr_hg_ppc_peak_vals = [fsi_near_lfr_hg_ppc_peak_vals, near_lfr_hg_ppc_pk_val];
                    fsi_near_lfr_lg_ppc_peak_vals = [fsi_near_lfr_lg_ppc_peak_vals, near_lfr_lg_ppc_pk_val];
                    fsi_near_hfr_hg_ppc_peak_vals = [fsi_near_hfr_hg_ppc_peak_vals, near_hfr_hg_ppc_pk_val];
                    fsi_near_hfr_lg_ppc_peak_vals = [fsi_near_hfr_lg_ppc_peak_vals, near_hfr_lg_ppc_pk_val];
                    fsi_near_p1_hg_sts_peak_vals = [fsi_near_p1_hg_sts_peak_vals, near_p1_hg_sts_pk_val];
                    fsi_near_p1_lg_sts_peak_vals = [fsi_near_p1_lg_sts_peak_vals, near_p1_lg_sts_pk_val];
                    fsi_near_p2_hg_sts_peak_vals = [fsi_near_p2_hg_sts_peak_vals, near_p2_hg_sts_pk_val];
                    fsi_near_p2_lg_sts_peak_vals = [fsi_near_p2_lg_sts_peak_vals, near_p2_lg_sts_pk_val];
                    fsi_near_p1_hg_ppc_peak_vals = [fsi_near_p1_hg_ppc_peak_vals, near_p1_hg_ppc_pk_val];
                    fsi_near_p1_lg_ppc_peak_vals = [fsi_near_p1_lg_ppc_peak_vals, near_p1_lg_ppc_pk_val];
                    fsi_near_p2_hg_ppc_peak_vals = [fsi_near_p2_hg_ppc_peak_vals, near_p2_hg_ppc_pk_val];
                    fsi_near_p2_lg_ppc_peak_vals = [fsi_near_p2_lg_ppc_peak_vals, near_p2_lg_ppc_pk_val];
                    
                    fsi_near_lfr_sts_corrs = [fsi_near_lfr_sts_corrs, near_lfr_sts_corr];
                    fsi_near_hfr_sts_corrs = [fsi_near_hfr_sts_corrs, near_hfr_sts_corr];
                    fsi_near_lfr_ppc_corrs = [fsi_near_lfr_ppc_corrs, near_lfr_ppc_corr];
                    fsi_near_hfr_ppc_corrs = [fsi_near_hfr_ppc_corrs, near_hfr_ppc_corr];
                    fsi_near_frs = [fsi_near_frs, near_spec_mean];
                    fsi_near_lfr_frs = [fsi_near_lfr_frs, near_lfr_spec_mean];
                    fsi_near_hfr_frs = [fsi_near_hfr_frs, near_hfr_spec_mean];
                    fsi_near_p1_sts_corrs = [fsi_near_p1_sts_corrs, near_p1_sts_corr];
                    fsi_near_p2_sts_corrs = [fsi_near_p2_sts_corrs, near_p2_sts_corr];
                    fsi_near_p1_ppc_corrs = [fsi_near_p1_ppc_corrs, near_p1_ppc_corr];
                    fsi_near_p2_ppc_corrs = [fsi_near_p2_ppc_corrs, near_p2_ppc_corr];
                    fsi_near_p1_frs = [fsi_near_p1_frs, near_p1_spec_mean];
                    fsi_near_p2_frs = [fsi_near_p2_frs, near_p2_spec_mean];
                end
            else
                continue
            end 
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        msn_near_sts_lg_peak_vals = [];
        msn_near_sts_hg_peak_vals = [];
        msn_near_lfr_hg_sts_peak_vals = [];
        msn_near_lfr_lg_sts_peak_vals = [];
        msn_near_hfr_hg_sts_peak_vals = [];
        msn_near_hfr_lg_sts_peak_vals = [];
        msn_near_ppc_lg_peak_vals = [];
        msn_near_ppc_hg_peak_vals = [];
        msn_near_lfr_hg_ppc_peak_vals = [];
        msn_near_lfr_lg_ppc_peak_vals = [];
        msn_near_hfr_hg_ppc_peak_vals = [];
        msn_near_hfr_lg_ppc_peak_vals = [];
        msn_near_p1_hg_sts_peak_vals = [];
        msn_near_p1_lg_sts_peak_vals = [];
        msn_near_p2_hg_sts_peak_vals = [];
        msn_near_p2_lg_sts_peak_vals = [];
        msn_near_p1_hg_ppc_peak_vals = [];
        msn_near_p1_lg_ppc_peak_vals = [];
        msn_near_p2_hg_ppc_peak_vals = [];
        msn_near_p2_lg_ppc_peak_vals = [];
        
        for iC = 1:length(msn_labels)
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
                pks = findpeaks(od.msn_res.near_p1_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_p1_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p1_lg_sts_peak = false;
                    end
                else
                    flag_near_p1_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_p1_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_p1_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p1_lg_ppc_peak = false;
                    end
                else
                    flag_near_p1_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_p1_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_p1_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_p1_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_p1_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p1_hg_sts_peak = false;
                    end
                else
                    flag_near_p1_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_p1_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_p1_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p1_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p1_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p1_hg_ppc_peak = false;
                    end
                else
                    flag_near_p1_hg_ppc_peak = false; 
                end
                
                % Peaks in near_p2_spec
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_p2_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_p2_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_p2_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_p2_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_lg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p2_lg_sts_peak = false;
                    end
                else
                    flag_near_p2_lg_sts_peak = false; 
                end                
                % Find Low Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_p2_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_p2_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_lg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p2_lg_ppc_peak = false;
                    end
                else
                    flag_near_p2_lg_ppc_peak = false; 
                end    
                % Find High Gamma STS Peak
                lf = find(od.msn_res.near_p2_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_p2_spec{iC}.freqs <= hg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_p2_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_p2_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_hg_sts_pk = pks.loc(loc1);
                    else
                       flag_near_p2_hg_sts_peak = false;
                    end
                else
                    flag_near_p2_hg_sts_peak = false; 
                end    
                % Find High Gamma PPC Peak
                pks = findpeaks(od.msn_res.near_p2_spec{iC}.ppc(lf:rf));
                max_val = max(od.msn_res.near_p2_spec{iC}.ppc);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_p2_spec{iC}.ppc(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        near_p2_hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_p2_hg_ppc_peak = false;
                    end
                else
                    flag_near_p2_hg_ppc_peak = false; 
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
                near_lfr_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,od.msn_res.near_lfr_spec{iC}.sts_vals);
                near_lfr_sts_corr = near_lfr_sts_corr(1,2);
                near_hfr_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,od.msn_res.near_hfr_spec{iC}.sts_vals);
                near_hfr_sts_corr = near_hfr_sts_corr(1,2);
                near_lfr_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,od.msn_res.near_lfr_spec{iC}.ppc);
                near_lfr_ppc_corr = near_lfr_ppc_corr(1,2);
                near_hfr_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,od.msn_res.near_hfr_spec{iC}.ppc);
                near_hfr_ppc_corr = near_hfr_ppc_corr(1,2);
                 
                p1_trials = od.msn_res.near_spec{iC}.p1_trials;     
                p2_trials = od.msn_res.near_spec{iC}.p2_trials;  
                p1_trials = p1_trials & nz_trials';
                p2_trials = p2_trials & nz_trials';
                near_p1_spec_mean = sum(od.msn_res.near_spec{iC}.trial_spk_count(p1_trials))/sum(od.msn_res.near_spec{iC}.trial_spk_count(p1_trials)'./od.msn_res.near_spec{iC}.mfr(p1_trials));
                near_p1_spec_range = [min(od.msn_res.near_spec{iC}.mfr(p1_trials)), max(od.msn_res.near_spec{iC}.mfr(p1_trials))];
                near_p2_spec_mean = sum(od.msn_res.near_spec{iC}.trial_spk_count(p2_trials))/sum(od.msn_res.near_spec{iC}.trial_spk_count(p2_trials)'./od.msn_res.near_spec{iC}.mfr(p2_trials));
                near_p2_spec_range = [min(od.msn_res.near_spec{iC}.mfr(p2_trials)), max(od.msn_res.near_spec{iC}.mfr(p2_trials))];
                near_p1_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,od.msn_res.near_p1_spec{iC}.sts_vals);
                near_p1_sts_corr = near_p1_sts_corr(1,2);
                near_p2_sts_corr  = corrcoef(od.msn_res.near_spec{iC}.sts_vals ,od.msn_res.near_p2_spec{iC}.sts_vals);
                near_p2_sts_corr = near_p2_sts_corr(1,2);
                near_p1_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,od.msn_res.near_p1_spec{iC}.ppc);
                near_p1_ppc_corr = near_p1_ppc_corr(1,2);
                near_p2_ppc_corr  = corrcoef(od.msn_res.near_spec{iC}.ppc ,od.msn_res.near_p2_spec{iC}.ppc);
                near_p2_ppc_corr = near_p2_ppc_corr(1,2);
                
                % if peaks exist in both low gamma and high gamma range
                if all_sts_peaks_in_lg && all_ppc_peaks_in_lg && ...
                    all_sts_peaks_in_hg && all_ppc_peaks_in_hg ...
                    
                    near_lg_sts_pk_val = od.msn_res.near_spec{iC}.sts_vals(lg(1)+near_lg_sts_pk-1);
                    near_hg_sts_pk_val = od.msn_res.near_spec{iC}.sts_vals(hg(1)+near_hg_sts_pk-1);
                    near_lfr_lg_sts_pk_val = od.msn_res.near_lfr_spec{iC}.sts_vals(lg(1)+near_lfr_lg_sts_pk-1);
                    near_lfr_hg_sts_pk_val = od.msn_res.near_lfr_spec{iC}.sts_vals(hg(1)+near_lfr_hg_sts_pk-1);
                    near_hfr_lg_sts_pk_val = od.msn_res.near_hfr_spec{iC}.sts_vals(lg(1)+near_hfr_lg_sts_pk-1);
                    near_hfr_hg_sts_pk_val = od.msn_res.near_hfr_spec{iC}.sts_vals(hg(1)+near_hfr_hg_sts_pk-1);
                    near_lg_ppc_pk_val = od.msn_res.near_spec{iC}.ppc(lg(1)+near_lg_ppc_pk-1);
                    near_hg_ppc_pk_val = od.msn_res.near_spec{iC}.ppc(hg(1)+near_hg_ppc_pk-1);
                    near_lfr_lg_ppc_pk_val = od.msn_res.near_lfr_spec{iC}.ppc(lg(1)+near_lfr_lg_ppc_pk-1);
                    near_lfr_hg_ppc_pk_val = od.msn_res.near_lfr_spec{iC}.ppc(hg(1)+near_lfr_hg_ppc_pk-1);
                    near_hfr_lg_ppc_pk_val = od.msn_res.near_hfr_spec{iC}.ppc(lg(1)+near_hfr_lg_ppc_pk-1);
                    near_hfr_hg_ppc_pk_val = od.msn_res.near_hfr_spec{iC}.ppc(hg(1)+near_hfr_hg_ppc_pk-1);  
                    near_p1_lg_sts_pk_val = od.msn_res.near_p1_spec{iC}.sts_vals(lg(1)+near_p1_lg_sts_pk-1);
                    near_p1_hg_sts_pk_val = od.msn_res.near_p1_spec{iC}.sts_vals(hg(1)+near_p1_hg_sts_pk-1);
                    near_p2_lg_sts_pk_val = od.msn_res.near_p2_spec{iC}.sts_vals(lg(1)+near_p2_lg_sts_pk-1);
                    near_p2_hg_sts_pk_val = od.msn_res.near_p2_spec{iC}.sts_vals(hg(1)+near_p2_hg_sts_pk-1);
                    near_p1_lg_ppc_pk_val = od.msn_res.near_p1_spec{iC}.ppc(lg(1)+near_p1_lg_ppc_pk-1);
                    near_p1_hg_ppc_pk_val = od.msn_res.near_p1_spec{iC}.ppc(hg(1)+near_p1_hg_ppc_pk-1);
                    near_p2_lg_ppc_pk_val = od.msn_res.near_p2_spec{iC}.ppc(lg(1)+near_p2_lg_ppc_pk-1);
                    near_p2_hg_ppc_pk_val = od.msn_res.near_p2_spec{iC}.ppc(hg(1)+near_p2_hg_ppc_pk-1);
                    
                    near_lg_sts_pk = od.msn_res.near_spec{iC}.freqs(lg(1)+near_lg_sts_pk-1);
                    near_hg_sts_pk = od.msn_res.near_spec{iC}.freqs(hg(1)+near_hg_sts_pk-1);
                    near_lfr_lg_sts_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_sts_pk-1);
                    near_lfr_hg_sts_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_sts_pk-1);
                    near_hfr_lg_sts_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_sts_pk-1);
                    near_hfr_hg_sts_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_sts_pk-1);
                    near_lg_ppc_pk = od.msn_res.near_spec{iC}.freqs(lg(1)+near_lg_ppc_pk-1);
                    near_hg_ppc_pk = od.msn_res.near_spec{iC}.freqs(hg(1)+near_hg_ppc_pk-1);
                    near_lfr_lg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
                    near_lfr_hg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
                    near_hfr_lg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
                    near_hfr_hg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);  
                    near_p1_lg_sts_pk = od.msn_res.near_p1_spec{iC}.freqs(lg(1)+near_p1_lg_sts_pk-1);
                    near_p1_hg_sts_pk = od.msn_res.near_p1_spec{iC}.freqs(hg(1)+near_p1_hg_sts_pk-1);
                    near_p2_lg_sts_pk = od.msn_res.near_p2_spec{iC}.freqs(lg(1)+near_p2_lg_sts_pk-1);
                    near_p2_hg_sts_pk = od.msn_res.near_p2_spec{iC}.freqs(hg(1)+near_p2_hg_sts_pk-1);
                    near_p1_lg_ppc_pk = od.msn_res.near_p1_spec{iC}.freqs(lg(1)+near_p1_lg_ppc_pk-1);
                    near_p1_hg_ppc_pk = od.msn_res.near_p1_spec{iC}.freqs(hg(1)+near_p1_hg_ppc_pk-1);
                    near_p2_lg_ppc_pk = od.msn_res.near_p2_spec{iC}.freqs(lg(1)+near_p2_lg_ppc_pk-1);
                    near_p2_hg_ppc_pk = od.msn_res.near_p2_spec{iC}.freqs(hg(1)+near_p2_hg_ppc_pk-1);
                    
                    msn_near_sts_lg_peak_vals = [msn_near_sts_lg_peak_vals, near_lg_sts_pk_val];
                    msn_near_sts_hg_peak_vals = [msn_near_sts_hg_peak_vals, near_hg_sts_pk_val];
                    msn_near_lfr_hg_sts_peak_vals = [msn_near_lfr_hg_sts_peak_vals, near_lfr_hg_sts_pk_val];
                    msn_near_lfr_lg_sts_peak_vals = [msn_near_lfr_lg_sts_peak_vals, near_lfr_lg_sts_pk_val];
                    msn_near_hfr_hg_sts_peak_vals = [msn_near_hfr_hg_sts_peak_vals, near_hfr_hg_sts_pk_val];
                    msn_near_hfr_lg_sts_peak_vals = [msn_near_hfr_lg_sts_peak_vals, near_hfr_lg_sts_pk_val];
                    msn_near_ppc_lg_peak_vals = [msn_near_ppc_lg_peak_vals, near_lg_ppc_pk_val];
                    msn_near_ppc_hg_peak_vals = [msn_near_ppc_hg_peak_vals, near_hg_ppc_pk_val];
                    msn_near_lfr_hg_ppc_peak_vals = [msn_near_lfr_hg_ppc_peak_vals, near_lfr_hg_ppc_pk_val];
                    msn_near_lfr_lg_ppc_peak_vals = [msn_near_lfr_lg_ppc_peak_vals, near_lfr_lg_ppc_pk_val];
                    msn_near_hfr_hg_ppc_peak_vals = [msn_near_hfr_hg_ppc_peak_vals, near_hfr_hg_ppc_pk_val];
                    msn_near_hfr_lg_ppc_peak_vals = [msn_near_hfr_lg_ppc_peak_vals, near_hfr_lg_ppc_pk_val];
                    msn_near_p1_hg_sts_peak_vals = [msn_near_p1_hg_sts_peak_vals, near_p1_hg_sts_pk_val];
                    msn_near_p1_lg_sts_peak_vals = [msn_near_p1_lg_sts_peak_vals, near_p1_lg_sts_pk_val];
                    msn_near_p2_hg_sts_peak_vals = [msn_near_p2_hg_sts_peak_vals, near_p2_hg_sts_pk_val];
                    msn_near_p2_lg_sts_peak_vals = [msn_near_p2_lg_sts_peak_vals, near_p2_lg_sts_pk_val];
                    msn_near_p1_hg_ppc_peak_vals = [msn_near_p1_hg_ppc_peak_vals, near_p1_hg_ppc_pk_val];
                    msn_near_p1_lg_ppc_peak_vals = [msn_near_p1_lg_ppc_peak_vals, near_p1_lg_ppc_pk_val];
                    msn_near_p2_hg_ppc_peak_vals = [msn_near_p2_hg_ppc_peak_vals, near_p2_hg_ppc_pk_val];
                    msn_near_p2_lg_ppc_peak_vals = [msn_near_p2_lg_ppc_peak_vals, near_p2_lg_ppc_pk_val];
                    
                    msn_near_sts_lg_peaks = [msn_near_sts_lg_peaks, near_lg_sts_pk];
                    msn_near_sts_hg_peaks = [msn_near_sts_hg_peaks, near_hg_sts_pk];
                    msn_near_lfr_hg_sts_peaks = [msn_near_lfr_hg_sts_peaks, near_lfr_hg_sts_pk];
                    msn_near_lfr_lg_sts_peaks = [msn_near_lfr_lg_sts_peaks, near_lfr_lg_sts_pk];
                    msn_near_hfr_hg_sts_peaks = [msn_near_hfr_hg_sts_peaks, near_hfr_hg_sts_pk];
                    msn_near_hfr_lg_sts_peaks = [msn_near_hfr_lg_sts_peaks, near_hfr_lg_sts_pk];
                    msn_near_ppc_lg_peaks = [msn_near_ppc_lg_peaks, near_lg_ppc_pk];
                    msn_near_ppc_hg_peaks = [msn_near_ppc_hg_peaks, near_hg_ppc_pk];
                    msn_near_lfr_hg_ppc_peaks = [msn_near_lfr_hg_ppc_peaks, near_lfr_hg_ppc_pk];
                    msn_near_lfr_lg_ppc_peaks = [msn_near_lfr_lg_ppc_peaks, near_lfr_lg_ppc_pk];
                    msn_near_hfr_hg_ppc_peaks = [msn_near_hfr_hg_ppc_peaks, near_hfr_hg_ppc_pk];
                    msn_near_hfr_lg_ppc_peaks = [msn_near_hfr_lg_ppc_peaks, near_hfr_lg_ppc_pk]; 
                    msn_near_p1_hg_sts_peaks = [msn_near_p1_hg_sts_peaks, near_p1_hg_sts_pk];
                    msn_near_p1_lg_sts_peaks = [msn_near_p1_lg_sts_peaks, near_p1_lg_sts_pk];
                    msn_near_p2_hg_sts_peaks = [msn_near_p2_hg_sts_peaks, near_p2_hg_sts_pk];
                    msn_near_p2_lg_sts_peaks = [msn_near_p2_lg_sts_peaks, near_p2_lg_sts_pk];
                    msn_near_p1_hg_ppc_peaks = [msn_near_p1_hg_ppc_peaks, near_p1_hg_ppc_pk];
                    msn_near_p1_lg_ppc_peaks = [msn_near_p1_lg_ppc_peaks, near_p1_lg_ppc_pk];
                    msn_near_p2_hg_ppc_peaks = [msn_near_p2_hg_ppc_peaks, near_p2_hg_ppc_pk];
                    msn_near_p2_lg_ppc_peaks = [msn_near_p2_lg_ppc_peaks, near_p2_lg_ppc_pk];
                    
                    msn_near_lfr_sts_corrs = [msn_near_lfr_sts_corrs, near_lfr_sts_corr];
                    msn_near_hfr_sts_corrs = [msn_near_hfr_sts_corrs, near_hfr_sts_corr];
                    msn_near_lfr_ppc_corrs = [msn_near_lfr_ppc_corrs, near_lfr_ppc_corr];
                    msn_near_hfr_ppc_corrs = [msn_near_hfr_ppc_corrs, near_hfr_ppc_corr];
                    msn_near_frs = [msn_near_frs, near_spec_mean];
                    msn_near_lfr_frs = [msn_near_lfr_frs, near_lfr_spec_mean];
                    msn_near_hfr_frs = [msn_near_hfr_frs, near_hfr_spec_mean];
                    msn_near_p1_sts_corrs = [msn_near_p1_sts_corrs, near_p1_sts_corr];
                    msn_near_p2_sts_corrs = [msn_near_p2_sts_corrs, near_p2_sts_corr];
                    msn_near_p1_ppc_corrs = [msn_near_p1_ppc_corrs, near_p1_ppc_corr];
                    msn_near_p2_ppc_corrs = [msn_near_p2_ppc_corrs, near_p2_ppc_corr];
                    msn_near_p1_frs = [msn_near_p1_frs, near_p1_spec_mean];
                    msn_near_p2_frs = [msn_near_p2_frs, near_p2_spec_mean];
                end
            else
                continue
            end
        end
        dummy = 1;
        if max(msn_near_ppc_lg_peak_vals) >= min(fsi_near_ppc_lg_peak_vals) | ...
           max(msn_near_hfr_lg_ppc_peak_vals) >= min(fsi_near_hfr_lg_ppc_peak_vals) | ...
           max(msn_near_lfr_lg_ppc_peak_vals) >= min(fsi_near_lfr_lg_ppc_peak_vals) | ...
           max(msn_near_ppc_hg_peak_vals) >= min(fsi_near_ppc_hg_peak_vals) | ...
           max(msn_near_hfr_hg_ppc_peak_vals) >= min(fsi_near_hfr_hg_ppc_peak_vals) | ...
           max(msn_near_lfr_hg_ppc_peak_vals) >= min(fsi_near_lfr_hg_ppc_peak_vals)
            
                sessions_with_msn_ppc_peak_ge_fsi = [sessions_with_msn_ppc_peak_ge_fsi; extractBefore(ofiles(jdx).name, '_')];
        end
    end
end

%% Plot Distance betweem Peak Frequencies vs Difference in Mean Firing rates
figure;
subplot(2,2,1);
scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_near_hfr_lg_ppc_peaks - fsi_near_lfr_lg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_near_hfr_hg_ppc_peaks - fsi_near_lfr_hg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('FSI PPC Peak Dif')
xlabel('FSI Firing Rate Dif')

subplot(2,2,3);
scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_lg_ppc_peaks - msn_near_lfr_lg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_hg_ppc_peaks - msn_near_lfr_hg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('MSN PPC Peak Dif')
xlabel('MSN Firing Rate Dif')

subplot(2,2,2);
scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_near_hfr_lg_sts_peaks - fsi_near_lfr_lg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_near_hfr_hg_sts_peaks - fsi_near_lfr_hg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('FSI STS Peak Dif')
xlabel('FSI Firing Rate Dif')

subplot(2,2,4);
scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_lg_sts_peaks - msn_near_lfr_lg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_hg_sts_peaks - msn_near_lfr_hg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('MSN STS Peak Dif')
xlabel('MSN Firing Rate Dif')


%% Plot Distance betweem Peak Frequencies vs Difference in Mean Firing rates for the control split
figure;
subplot(2,2,1);
scatter(fsi_near_p2_frs - fsi_near_p1_frs, fsi_near_p2_lg_ppc_peaks - fsi_near_p1_lg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(fsi_near_p2_frs - fsi_near_p1_frs, fsi_near_p2_hg_ppc_peaks - fsi_near_p1_hg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('FSI PPC Peak Dif')
xlabel('FSI Firing Rate Dif')

subplot(2,2,3);
scatter(msn_near_p2_frs - msn_near_p1_frs, msn_near_p2_lg_ppc_peaks - msn_near_p1_lg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(msn_near_p2_frs - msn_near_p1_frs, msn_near_p2_hg_ppc_peaks - msn_near_p1_hg_ppc_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('MSN PPC Peak Dif')
xlabel('MSN Firing Rate Dif')

subplot(2,2,2);
scatter(fsi_near_p2_frs - fsi_near_p1_frs, fsi_near_p2_lg_sts_peaks - fsi_near_p1_lg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(fsi_near_p2_frs - fsi_near_p1_frs, fsi_near_p2_hg_sts_peaks - fsi_near_p1_hg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('FSI STS Peak Dif')
xlabel('FSI Firing Rate Dif')

subplot(2,2,4);
scatter(msn_near_p2_frs - msn_near_p1_frs, msn_near_p2_lg_sts_peaks - msn_near_p1_lg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
hold on;
scatter(msn_near_p2_frs - msn_near_p1_frs, msn_near_p2_hg_sts_peaks - msn_near_p1_hg_sts_peaks, ...
    'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
ylabel('MSN STS Peak Dif')
xlabel('MSN Firing Rate Dif')

%% Plot distribution of correlations
figure
subplot(2,2,1)
hold on
histogram(fsi_near_hfr_ppc_corrs, -0.05:0.1:1.05,'FaceColor','green','FaceAlpha',0.4);
histogram(fsi_near_lfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('FSI PPC Correlation');
subplot(2,2,2)
hold on
histogram(fsi_near_hfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(fsi_near_lfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('FSI STS Correlation');
subplot(2,2,3)
hold on
histogram(msn_near_hfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(msn_near_lfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('MSN PPC Correlation');
subplot(2,2,4)
hold on
histogram(msn_near_hfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(msn_near_lfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('MSN STS Correlation');

%% Plot distribution of correlations
figure
subplot(2,2,1)
hold on
histogram(fsi_near_hfr_ppc_corrs, -0.05:0.1:1.05,'FaceColor','green','FaceAlpha',0.4);
histogram(fsi_near_lfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('FSI PPC Correlation');
subplot(2,2,2)
hold on
histogram(fsi_near_hfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(fsi_near_lfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('FSI STS Correlation');
subplot(2,2,3)
hold on
histogram(msn_near_hfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(msn_near_lfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('MSN PPC Correlation');
subplot(2,2,4)
hold on
histogram(msn_near_hfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(msn_near_lfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('MSN STS Correlation');

%% Plot distribution of correlations of the control split
figure
subplot(2,2,1)
hold on
histogram(fsi_near_hfr_ppc_corrs, -0.05:0.1:1.05,'FaceColor','green','FaceAlpha',0.4);
histogram(fsi_near_lfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('FSI PPC Correlation');
subplot(2,2,2)
hold on
histogram(fsi_near_hfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(fsi_near_lfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('FSI STS Correlation');
subplot(2,2,3)
hold on
histogram(msn_near_hfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(msn_near_lfr_ppc_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('MSN PPC Correlation');
subplot(2,2,4)
hold on
histogram(msn_near_hfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','green','FaceAlpha',0.4);
histogram(msn_near_lfr_sts_corrs, -0.05:0.1:1.05, 'FaceColor','red','FaceAlpha',0.4)
title('MSN STS Correlation');

%% Plot Scatter plots of difference between PPC peak frequncies vs difference in mean firing rate, color coded by correlation
corr_thresh = 0.65;
mask_hfr_corrs = (msn_near_hfr_ppc_corrs < corr_thresh);
mask_lfr_corrs = (msn_near_lfr_ppc_corrs < corr_thresh);
figure;
subplot(2,1,1)
scatter(msn_near_hfr_frs(mask_hfr_corrs & mask_lfr_corrs) - msn_near_lfr_frs(mask_hfr_corrs & mask_lfr_corrs), msn_near_hfr_lg_ppc_peaks(mask_hfr_corrs & mask_lfr_corrs) - msn_near_lfr_lg_ppc_peaks(mask_hfr_corrs & mask_lfr_corrs), ...
'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
% subplot(2,1,2)
hold on
scatter(msn_near_hfr_frs(~mask_hfr_corrs & ~mask_lfr_corrs) - msn_near_lfr_frs(~mask_hfr_corrs & ~mask_lfr_corrs), msn_near_hfr_lg_ppc_peaks(~mask_hfr_corrs & ~mask_lfr_corrs) - msn_near_lfr_lg_ppc_peaks(~mask_hfr_corrs & ~mask_lfr_corrs), ...
'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);