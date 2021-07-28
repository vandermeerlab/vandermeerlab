cd('D:\RandomVstrAnalysis\ft_results');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
lg = [30,65];
hg = [65,100];
pk_thresh = 0.5;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        
        msn_labels  = od.label(od.cell_type == 1);
        for iC = 1:length(msn_labels)
            o_prefix = extractBefore(msn_labels{iC},'.t');
            
            % Plot On Track stuff first
            fig = figure('WindowState', 'maximized');
            subplot(2,3,1)
            plot(od.msn_res.onTrack_spec{iC}.sta_time, od.msn_res.onTrack_spec{iC}.sta_vals);
            title(sprintf('On Track STA, %d spikes', od.msn_res.onTrack_spec{iC}.spk_count));
            
            if ~od.msn_res.onTrack_spec{iC}.flag_nansts
                subplot(2,3,2)
                plot(od.msn_res.onTrack_spec{iC}.freqs, od.msn_res.onTrack_spec{iC}.sts_vals);
                xlabel('Freqs')
                title('On Track STS');
            end
            
            if ~od.msn_res.onTrack_spec{iC}.flag_nanppc
                subplot(2,3,3)
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
                subplot(2,3,4)
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                flag_leg = true;
                this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.msn_res.near_spec{iC}.spk_count);
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.msn_res.near_lfr_spec{iC}.spk_count);
                plot(od.msn_res.near_spec{iC}.sta_time, od.msn_res.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.msn_res.near_hfr_spec{iC}.spk_count);
                legend(this_legend, 'Location', 'southwest')
                title('Near Reward STA');
                
                flag_near_lg_sts_peak = true;
                flag_near_hg_sts_peak = true;
                flag_near_lg_ppc_peak = true;
                flag_near_hg_ppc_peak = true;
                
                % Find Low Gamma STS Peak
                lf = find(od.msn_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                rf = find(od.msn_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
                pks = findpeaks(od.msn_res.near_spec{iC}.sts_vals(lf:rf));
                max_val = max(od.msn_res.near_spec{iC}.sts_vals);
                if ~isempty(pks.loc) 
                    [pk1, loc1] = max(od.msn_res.near_spec{iC}.sts_vals(lf+pks.loc-1));
                    if pk1 > max_val*pk_thresh
                        lg_sts_pk = pks.loc(loc1);
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
                        lg_ppc_pk = pks.loc(loc1);
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
                        hg_sts_pk = pks.loc(loc1);
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
                        hg_ppc_pk = pks.loc(loc1);
                    else
                       flag_near_hg_ppc_peak = false;
                    end
                else
                    flag_near_hg_ppc_peak = false; 
                end
                
                % Plot STS_stuff first    
                if flag_lg_ppc_peak 
                
                flag_leg = false;
                subplot(2,3,5)
                if ~od.msn_res.near_spec{iC}.flag_nansts && od.msn_res.near_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.msn_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.msn_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.msn_res.near_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_spec{iC}.sts_vals(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_sts_lg_peak = od.msn_res.near_spec{iC}.freqs(lf+p1-1);
                    plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.sts_vals, 'blue');
                    q0 = xline(near_sts_lg_peak, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    hold on;
                    % Plot high gamma peak
                    lf = find(od.msn_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.msn_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.msn_res.near_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_spec{iC}.sts_vals(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_sts_hg_peak = od.msn_res.near_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_sts_hg_peak, '--blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';      
                    this_legend = {};
                    flag_leg = true;
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_sts_lg_peak, near_sts_hg_peak);
                else
                    close all;
                    continue;
                end
                if ~od.msn_res.near_lfr_spec{iC}.flag_nansts && od.msn_res.near_lfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.msn_res.near_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.msn_res.near_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_lfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.msn_res.near_lfr_spec{iC}.sts_vals);
                     if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_lfr_sts_lg_peak = od.msn_res.near_lfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.sts_vals, 'Color', 'red');
                    q0 = xline(near_lfr_sts_lg_peak, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.msn_res.near_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.msn_res.near_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_lfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.msn_res.near_lfr_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_lfr_sts_hg_peak = od.msn_res.near_lfr_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_lfr_sts_hg_peak, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('LFR trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_lfr_sts_lg_peak, near_lfr_sts_hg_peak);   
                else
                    close all;
                    continue;
                end
                if ~od.msn_res.near_hfr_spec{iC}.flag_nansts && od.msn_res.near_hfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.msn_res.near_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.msn_res.near_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_hfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.msn_res.near_hfr_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_hfr_sts_lg_peak = od.msn_res.near_hfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.sts_vals, 'Color', 'green');
                    q0 = xline(near_hfr_sts_lg_peak, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.msn_res.near_hfr_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.msn_res.near_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_hfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.msn_res.near_hfr_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_hfr_sts_hg_peak = od.msn_res.near_hfr_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_hfr_sts_hg_peak, '--green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('HFR trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_hfr_sts_lg_peak, near_hfr_sts_hg_peak); 
                else
                    close all;
                    continue;
                end
                if flag_leg
                    legend(this_legend, 'Location', 'northwest')
                    xlabel('Freqs')
                    title('Near Reward STS');
                end
                
                flag_leg = false;
                subplot(2,3,6)
                if ~od.msn_res.near_spec{iC}.flag_nanppc && od.msn_res.near_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.msn_res.near_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.msn_res.near_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_spec{iC}.ppc(lf:rf));
                    max_val = max(od.msn_res.near_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_spec{iC}.ppc(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_ppc_lg_peak = od.msn_res.near_spec{iC}.freqs(lf+p1-1);
                    plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.ppc, 'blue');
                    q0 = xline(near_ppc_lg_peak, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    hold on;
                    % Plot high gamma peak
                    lf = find(od.msn_res.near_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.msn_res.near_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_spec{iC}.ppc(lf:rf));
                    max_val = max(od.msn_res.near_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_spec{iC}.ppc(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_ppc_hg_peak = od.msn_res.near_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_ppc_hg_peak, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend = {};
                    flag_leg = true;
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_ppc_lg_peak, near_ppc_hg_peak); 
                else
                    close all;
                    continue;
                end
                if ~od.msn_res.near_lfr_spec{iC}.flag_nanppc && od.msn_res.near_lfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.msn_res.near_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.msn_res.near_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_lfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.msn_res.near_lfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_lfr_ppc_lg_peak = od.msn_res.near_lfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.msn_res.near_lfr_spec{iC}.freqs, od.msn_res.near_lfr_spec{iC}.ppc, 'Color', 'red');
                    q0 = xline(near_lfr_ppc_lg_peak, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.msn_res.near_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.msn_res.near_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_lfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.msn_res.near_lfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_lfr_ppc_hg_peak = od.msn_res.near_lfr_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_lfr_ppc_hg_peak, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_lfr_ppc_lg_peak, near_lfr_ppc_hg_peak); 
                else
                    close all;
                    continue;
                end
                if ~od.msn_res.near_hfr_spec{iC}.flag_nanppc && od.msn_res.near_hfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.msn_res.near_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.msn_res.near_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_hfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.msn_res.near_hfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_hfr_ppc_lg_peak = od.msn_res.near_hfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.msn_res.near_hfr_spec{iC}.freqs, od.msn_res.near_hfr_spec{iC}.ppc, 'Color', 'green');
                    q0 = xline(near_hfr_ppc_lg_peak, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.msn_res.near_hfr_spec{iC}.freqs >=hg(1), 1, 'first');
                    rf = find(od.msn_res.near_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.msn_res.near_hfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.msn_res.near_hfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.msn_res.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
                        if pk1 > max_val*pk_thresh
                            p1 = pks.loc(loc1);
                        else
                           close all;
                           continue;
                        end
                    else
                        close all;
                        continue;
                    end
                    near_hfr_ppc_hg_peak = od.msn_res.near_hfr_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_hfr_ppc_hg_peak, '--green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('HFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_hfr_ppc_lg_peak, near_hfr_ppc_hg_peak);
                else
                    close all;
                    continue;
                end
                if flag_leg
                    legend(this_legend, 'Location', 'northeast')
                    xlabel('Freqs')
                    title('Near Reward PPC');
                end
            else
                close all;
                continue;
            end
            

            
            % Indicate lg category of cell in suffix
            if near_hfr_ppc_lg_peak - near_lfr_ppc_lg_peak > 10
                o_prefix = cat(2, o_prefix, '_lg_ppc_g10');
            elseif near_hfr_ppc_lg_peak - near_lfr_ppc_lg_peak < -10
                o_prefix = cat(2, o_prefix, '_lg_ppc_gn10');
            else
                o_prefix = cat(2, o_prefix, '_lg_ppc_l10');
            end
            
            % Indicate hg category of cell in suffix
            if near_hfr_ppc_hg_peak - near_lfr_ppc_hg_peak > 10
                o_prefix = cat(2, o_prefix, '_hg_ppc_g10');
            elseif near_hfr_ppc_hg_peak - near_lfr_ppc_hg_peak < -10
                o_prefix = cat(2, o_prefix, '_hg_ppc_gn10');
            else
                o_prefix = cat(2, o_prefix, '_hg_ppc_l10');
            end
            
            %Indicate cell type
            if od.msn_res.cell_type(iC) == 1
                o_name = cat(2, o_prefix, '_MSN');
            else
                 o_name = cat(2, o_prefix, '_FSI');
            end
            
           WriteFig(fig,o_name,1);
           close all;  
        end
    end
end
