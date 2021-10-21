% cd('D:\RandomVstrAnalysis\ft_results');
cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
lg = [5,30];
hg = [30,100];
pk_thresh = -1;
num_control_splits = 100;

msn_counter  = {};
fsi_counter = {};
% Collect 
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        
        fsi_labels  = od.label(od.cell_type == 2);
        for iC = 1:length(fsi_labels)
            o_prefix = extractBefore(fsi_labels{iC},'.t');
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if isfield(od.fsi_res.near_spec{iC},'flag_no_control_split') && ...
                        ~od.fsi_res.near_spec{iC}.flag_no_control_split          

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
                    continue;
                end
            else
                continue
            end
            fsi_counter{length(fsi_counter)+1} = o_prefix;
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        for iC = 1:length(msn_labels)
            o_prefix = extractBefore(msn_labels{iC},'.t');
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if isfield(od.msn_res.near_spec{iC},'flag_no_control_split') && ...
                    ~od.msn_res.near_spec{iC}.flag_no_control_split          
                
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
                    continue;
                end 
            else
                continue
            end
            msn_counter{length(msn_counter)+1} = o_prefix;
        end
    end
end