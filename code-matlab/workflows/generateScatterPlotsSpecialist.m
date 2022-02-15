% Generates Scatter Plots of Difference in Peak Frequencies vs difference
% in mean firing rates

fsi_near_lfr_hg_ppc_peaks = [];
fsi_near_lfr_lg_ppc_peaks = [];
fsi_near_hfr_hg_ppc_peaks = [];
fsi_near_hfr_lg_ppc_peaks = [];
fsi_near_lfr_frs = [];
fsi_near_hfr_frs = [];
valid_fsi_labels = [];

msn_near_lfr_hg_ppc_peaks = [];
msn_near_lfr_lg_ppc_peaks = [];
msn_near_hfr_hg_ppc_peaks = [];
msn_near_hfr_lg_ppc_peaks = [];
msn_near_lfr_frs = [];
msn_near_hfr_frs = [];
valid_msn_labels = [];

sessions_with_msn_ppc_peak_ge_fsi = [];

% cd('D:\RandomVstrAnalysis\ft_results');
cd('D:\RandomVstrAnalysis\ft_results');

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
                    valid_fsi_labels = [valid_fsi_labels fsi_labels(iC)];
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
                    valid_msn_labels = [valid_msn_labels msn_labels(iC)];
                end
            else
                continue
            end
        end
    end
end

%% Plot Distance betweem Peak Frequencies vs Difference in Mean Firing rates in MSNs and FSIs
close all;
fig_1 = figure('WindowState', 'maximized');
s1 = scatter(msn_near_hfr_frs - msn_near_lfr_frs, msn_near_hfr_hg_ppc_peaks - msn_near_lfr_hg_ppc_peaks);
s1.Marker = 'o';
s1.MarkerFaceColor = [0.8500 0.3250 0.0980];
s1.MarkerFaceAlpha = 0.5;
s1.MarkerEdgeColor = [0.8500 0.3250 0.0980];
s1.MarkerEdgeAlpha = 0.5;
s1.SizeData = 150;
hold on;
s2 = scatter(fsi_near_hfr_frs - fsi_near_lfr_frs, fsi_near_hfr_hg_ppc_peaks - fsi_near_lfr_hg_ppc_peaks);
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
%%
fig_2 = figure('WindowState', 'maximized');
hold on;
for iF = 1:length(valid_fsi_labels)
    text(fsi_near_hfr_frs(iF) - fsi_near_lfr_frs(iF), ...
        fsi_near_hfr_hg_ppc_peaks(iF) - fsi_near_lfr_hg_ppc_peaks(iF), ...
        valid_fsi_labels(iF), 'Interpreter', 'none');
end
xlim([0,10])
ylim([-100 100])
%%
fig_3 = figure('WindowState', 'maximized');
hold on;
for iF = 1:length(valid_msn_labels)
    text(msn_near_hfr_frs(iF) - msn_near_lfr_frs(iF), ...
        msn_near_hfr_hg_ppc_peaks(iF) - msn_near_lfr_hg_ppc_peaks(iF), ...
        valid_msn_labels(iF), 'Interpreter', 'none');
end
xlim([0,10])
ylim([-100 100])
%%
WriteFig(fig,'D:\RandomVstrAnalysis\temp\test',1);
