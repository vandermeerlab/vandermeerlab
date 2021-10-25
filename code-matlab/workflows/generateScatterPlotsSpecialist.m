fsi_near_lfr_hg_ppc_peaks = [];
fsi_near_lfr_lg_ppc_peaks = [];
fsi_near_hfr_hg_ppc_peaks = [];
fsi_near_hfr_lg_ppc_peaks = [];
fsi_near_lfr_frs = [];
fsi_near_hfr_frs = [];

msn_near_lfr_hg_ppc_peaks = [];
msn_near_lfr_lg_ppc_peaks = [];
msn_near_hfr_hg_ppc_peaks = [];
msn_near_hfr_lg_ppc_peaks = [];
msn_near_lfr_frs = [];
msn_near_hfr_frs = [];


sessions_with_msn_ppc_peak_ge_fsi = [];

% cd('D:\RandomVstrAnalysis\ft_results');
cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
lg = [3,30];
hg = [5,100];
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
                
                flag_near_lfr_lg_ppc_peak = true;
                flag_near_lfr_hg_ppc_peak = true;
                flag_near_hfr_lg_ppc_peak = true;
                flag_near_hfr_hg_ppc_peak = true;

                % Peaks in near_lfr_spec
                % Find Low Gamma PPC Peak
                lf = find(od.fsi_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
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
                
                % Find High Gamma PPC Peak
                lf = find(od.fsi_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
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
                % Find Low Gamma PPC Peak
                lf = find(od.fsi_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
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
                
                % Find High Gamma PPC Peak
                lf = find(od.fsi_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.fsi_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
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
                
                all_ppc_peaks_in_lg = flag_near_hfr_lg_ppc_peak && ...
                    flag_near_lfr_lg_ppc_peak;
                         
                all_ppc_peaks_in_hg = flag_near_hfr_hg_ppc_peak && ...
                    flag_near_lfr_hg_ppc_peak;

                % Skip plotting only if no peaks in either hg or lg     
                if ~all_ppc_peaks_in_lg && ...
                      ~all_ppc_peaks_in_hg  
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
                
                % if peaks exist in both low gamma and high gamma range
                if all_ppc_peaks_in_lg && all_ppc_peaks_in_hg ...
                   
                    near_lfr_lg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
                    near_lfr_hg_ppc_pk = od.fsi_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
                    near_hfr_lg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
                    near_hfr_hg_ppc_pk = od.fsi_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);
                    
                    fsi_near_lfr_hg_ppc_peaks = [fsi_near_lfr_hg_ppc_peaks, near_lfr_hg_ppc_pk];
                    fsi_near_lfr_lg_ppc_peaks = [fsi_near_lfr_lg_ppc_peaks, near_lfr_lg_ppc_pk];
                    fsi_near_hfr_hg_ppc_peaks = [fsi_near_hfr_hg_ppc_peaks, near_hfr_hg_ppc_pk];
                    fsi_near_hfr_lg_ppc_peaks = [fsi_near_hfr_lg_ppc_peaks, near_hfr_lg_ppc_pk];
                    
                    fsi_near_lfr_frs = [fsi_near_lfr_frs, near_lfr_spec_mean];
                    fsi_near_hfr_frs = [fsi_near_hfr_frs, near_hfr_spec_mean];
                end
            else
                continue
            end 
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC},'flag_no_control_split') && ...
                    ~od.msn_res.near_spec{iC}.flag_no_control_split  
                
                flag_near_lfr_lg_ppc_peak = true;
                flag_near_lfr_hg_ppc_peak = true;
                flag_near_hfr_lg_ppc_peak = true;
                flag_near_hfr_hg_ppc_peak = true;
                
                % Peaks in near_lfr_spec              
                % Find Low Gamma PPC Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
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
                
                % Find High Gamma PPC Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
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
                % Find Low Gamma PPC Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
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
                
                % Find High Gamma PPC Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
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
                

                all_ppc_peaks_in_lg = flag_near_hfr_lg_ppc_peak && ...
                    flag_near_lfr_lg_ppc_peak;
                         
                all_ppc_peaks_in_hg = flag_near_hfr_hg_ppc_peak && ...
                    flag_near_lfr_hg_ppc_peak;

                % Skip plotting only if no peaks in either hg or lg     
                if ~all_ppc_peaks_in_lg && ~all_ppc_peaks_in_hg  
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
                
                % if peaks exist in both low gamma and high gamma range
                if all_ppc_peaks_in_lg && all_ppc_peaks_in_hg
                    
                    near_lfr_lg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(lg(1)+near_lfr_lg_ppc_pk-1);
                    near_lfr_hg_ppc_pk = od.msn_res.near_lfr_spec{iC}.freqs(hg(1)+near_lfr_hg_ppc_pk-1);
                    near_hfr_lg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(lg(1)+near_hfr_lg_ppc_pk-1);
                    near_hfr_hg_ppc_pk = od.msn_res.near_hfr_spec{iC}.freqs(hg(1)+near_hfr_hg_ppc_pk-1);  
                    

                    msn_near_lfr_hg_ppc_peaks = [msn_near_lfr_hg_ppc_peaks, near_lfr_hg_ppc_pk];
                    msn_near_lfr_lg_ppc_peaks = [msn_near_lfr_lg_ppc_peaks, near_lfr_lg_ppc_pk];
                    msn_near_hfr_hg_ppc_peaks = [msn_near_hfr_hg_ppc_peaks, near_hfr_hg_ppc_pk];
                    msn_near_hfr_lg_ppc_peaks = [msn_near_hfr_lg_ppc_peaks, near_hfr_lg_ppc_pk]; 
                    
                    msn_near_lfr_frs = [msn_near_lfr_frs, near_lfr_spec_mean];
                    msn_near_hfr_frs = [msn_near_hfr_frs, near_hfr_spec_mean];
                end
            else
                continue
            end
        end
    end
end

%% Plot Distance betweem Peak Frequencies vs Difference in Mean Firing rates
close all
fig = figure('WindowState', 'maximized');
scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_hg_ppc_peaks - msn_near_lfr_hg_ppc_peaks, ...
    'Marker', 'o', 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', 0.5);
ylabel('Difference in Frequency of Maximum Phase Locking', 'FontSize', 16)
xlabel('Difference in Firing Rates between HFR trials and LFR trials', 'FontSize', 16)
box off;
%%
WriteFig(fig,'~/test',1);
