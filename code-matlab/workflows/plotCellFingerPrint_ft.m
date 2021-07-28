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
        dummy = 2;
        num_cells = length(od.label);
        for iC = 1:num_cells
            o_prefix = extractBefore(od.label{iC},'.t');
            % Plot On Track stuff first
            fig = figure('WindowState', 'maximized');
            subplot(2,3,1)
            plot(od.onTrack_spec{iC}.sta_time, od.onTrack_spec{iC}.sta_vals);
            title(sprintf('On Track STA, %d spikes', od.onTrack_spec{iC}.spk_count));
            
            if ~od.onTrack_spec{iC}.flag_nansts
                subplot(2,3,2)
                plot(od.onTrack_spec{iC}.freqs, od.onTrack_spec{iC}.sts_vals);
                xlabel('Freqs')
                title('On Track STS');
            end
            
            if ~od.onTrack_spec{iC}.flag_nanppc
                subplot(2,3,3)
                plot(od.onTrack_spec{iC}.freqs, od.onTrack_spec{iC}.ppc);
                xlabel('Freqs')
                title('On Track PPC');
            end
            
            % Plot Near Reward stuff next
            flag_leg = false;
            if ~od.near_spec{iC}.flag_zeroSpikes && od.near_spec{iC}.spk_count >= 100
                subplot(2,3,4)
                plot(od.near_spec{iC}.sta_time, od.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                flag_leg = true;
                this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.near_spec{iC}.spk_count);
                if ~od.near_lfr_spec{iC}.flag_zeroSpikes && od.near_lfr_spec{iC}.spk_count >= 100
                    plot(od.near_spec{iC}.sta_time, od.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                    this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.near_lfr_spec{iC}.spk_count);
                end
                if ~od.near_hfr_spec{iC}.flag_zeroSpikes && od.near_hfr_spec{iC}.spk_count >= 100
                    plot(od.near_spec{iC}.sta_time, od.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                    this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.near_hfr_spec{iC}.spk_count);
                end
                if flag_leg
                    legend(this_legend, 'Location', 'southwest')
                    title('Near Reward STA');
                end
                
                flag_leg = false;
                subplot(2,3,5)
                if ~od.near_spec{iC}.flag_nansts && od.near_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.near_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.near_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.near_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.near_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_spec{iC}.sts_vals(lf+pks.loc-1));
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
                    near_sts_lg_peak = od.near_spec{iC}.freqs(lf+p1-1);
                    plot(od.near_spec{iC}.freqs, od.near_spec{iC}.sts_vals, 'blue');
                    q0 = xline(near_sts_lg_peak, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    hold on;
                    % Plot high gamma peak
                    lf = find(od.near_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.near_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.near_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.near_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_spec{iC}.sts_vals(lf+pks.loc-1));
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
                    near_sts_hg_peak = od.near_spec{iC}.freqs(lf+p1-1);
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
                if ~od.near_lfr_spec{iC}.flag_nansts && od.near_lfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.near_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.near_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.near_lfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.near_lfr_spec{iC}.sts_vals);
                     if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
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
                    near_lfr_sts_lg_peak = od.near_lfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.near_lfr_spec{iC}.freqs, od.near_lfr_spec{iC}.sts_vals, 'Color', 'red');
                    q0 = xline(near_lfr_sts_lg_peak, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.near_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.near_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.near_lfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.near_lfr_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
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
                    near_lfr_sts_hg_peak = od.near_lfr_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_lfr_sts_hg_peak, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('LFR trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_lfr_sts_lg_peak, near_lfr_sts_hg_peak);   
                else
                    close all;
                    continue;
                end
                if ~od.near_hfr_spec{iC}.flag_nansts && od.near_hfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.near_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.near_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.near_hfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.near_hfr_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
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
                    near_hfr_sts_lg_peak = od.near_hfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.near_hfr_spec{iC}.freqs, od.near_hfr_spec{iC}.sts_vals, 'Color', 'green');
                    q0 = xline(near_hfr_sts_lg_peak, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.near_hfr_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.near_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.near_hfr_spec{iC}.sts_vals(lf:rf));
                    max_val = max(od.near_hfr_spec{iC}.sts_vals);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
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
                    near_hfr_sts_hg_peak = od.near_hfr_spec{iC}.freqs(lf+p1-1);
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
                if ~od.near_spec{iC}.flag_nanppc && od.near_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.near_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.near_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.near_spec{iC}.ppc(lf:rf));
                    max_val = max(od.near_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_spec{iC}.ppc(lf+pks.loc-1));
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
                    near_ppc_lg_peak = od.near_spec{iC}.freqs(lf+p1-1);
                    plot(od.near_spec{iC}.freqs, od.near_spec{iC}.ppc, 'blue');
                    q0 = xline(near_ppc_lg_peak, 'blue');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    hold on;
                    % Plot high gamma peak
                    lf = find(od.near_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.near_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.near_spec{iC}.ppc(lf:rf));
                    max_val = max(od.near_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_spec{iC}.ppc(lf+pks.loc-1));
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
                    near_ppc_hg_peak = od.near_spec{iC}.freqs(lf+p1-1);
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
                if ~od.near_lfr_spec{iC}.flag_nanppc && od.near_lfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.near_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.near_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.near_lfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.near_lfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
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
                    near_lfr_ppc_lg_peak = od.near_lfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.near_lfr_spec{iC}.freqs, od.near_lfr_spec{iC}.ppc, 'Color', 'red');
                    q0 = xline(near_lfr_ppc_lg_peak, 'red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.near_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
                    rf = find(od.near_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.near_lfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.near_lfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_lfr_spec{iC}.ppc(lf+pks.loc-1));
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
                    near_lfr_ppc_hg_peak = od.near_lfr_spec{iC}.freqs(lf+p1-1);
                    q0 = xline(near_lfr_ppc_hg_peak, '--red');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    this_legend{length(this_legend)+1} = sprintf ...
                        ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
                            near_lfr_ppc_lg_peak, near_lfr_ppc_hg_peak); 
                else
                    close all;
                    continue;
                end
                if ~od.near_hfr_spec{iC}.flag_nanppc && od.near_hfr_spec{iC}.spk_count >= 100
                    % Plot low gamma peak
                    lf = find(od.near_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
                    rf = find(od.near_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
                    pks = findpeaks(od.near_hfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.near_hfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
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
                    near_hfr_ppc_lg_peak = od.near_hfr_spec{iC}.freqs(lf+p1-1);
                    plot(od.near_hfr_spec{iC}.freqs, od.near_hfr_spec{iC}.ppc, 'Color', 'green');
                    q0 = xline(near_hfr_ppc_lg_peak, 'green');
                    q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    % Plot high gamma peak
                    lf = find(od.near_hfr_spec{iC}.freqs >=hg(1), 1, 'first');
                    rf = find(od.near_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
                    pks = findpeaks(od.near_hfr_spec{iC}.ppc(lf:rf));
                    max_val = max(od.near_hfr_spec{iC}.ppc);
                    if ~isempty(pks.loc) 
                        [pk1, loc1] = max(od.near_hfr_spec{iC}.ppc(lf+pks.loc-1));
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
                    near_hfr_ppc_hg_peak = od.near_hfr_spec{iC}.freqs(lf+p1-1);
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
            
%             % Plot Away Reward stuff next
%             flag_leg = false;
%             if ~od.away_spec{iC}.flag_zeroSpikes && od.away_spec{iC}.spk_count >= 100
%                 subplot(3,3,7)
%                 plot(od.away_spec{iC}.sta_time, od.away_spec{iC}.sta_vals);
%                 hold on;
%                 this_legend = {};
%                 flag_leg = true;
%                 this_legend{length(this_legend)+1} = sprintf('All Trials: %d spikes',od.away_spec{iC}.spk_count);
%                 if ~od.away_lfr_spec{iC}.flag_zeroSpikes && od.away_lfr_spec{iC}.spk_count >= 100
%                     plot(od.away_spec{iC}.sta_time, od.away_lfr_spec{iC}.sta_vals, 'Color', 'red');
%                     this_legend{length(this_legend)+1} = sprintf('LFR Trials: %d spikes',od.away_lfr_spec{iC}.spk_count);
%                 end
%                 if ~od.away_hfr_spec{iC}.flag_zeroSpikes && od.away_hfr_spec{iC}.spk_count >= 100
%                     plot(od.away_spec{iC}.sta_time, od.away_hfr_spec{iC}.sta_vals, 'Color', 'green');
%                     this_legend{length(this_legend)+1} = sprintf('HFR Trials: %d spikes',od.away_hfr_spec{iC}.spk_count);
%                 end
%                 if flag_leg
%                     legend(this_legend, 'Location', 'southwest')
%                     title('Away Reward STA');
%                 end
%                 
%                 flag_leg = false;
%                 subplot(3,3,8)
%                 if ~od.away_spec{iC}.flag_nansts && od.away_spec{iC}.spk_count >= 100
%                     % Plot low gamma peak
%                     lf = find(od.away_spec{iC}.freqs >= lg(1), 1, 'first');
%                     rf = find(od.away_spec{iC}.freqs <= lg(2), 1, 'last');
%                     pks = findpeaks(od.away_spec{iC}.sts_vals(lf:rf));
%                     max_val = max(od.away_spec{iC}.sts_vals);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_spec{iC}.sts_vals(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_spec{iC}.sts_vals(lf:rf));
%                     end
%                     away_sts_lg_peak = od.away_spec{iC}.freqs(lf+p1-1);
%                     plot(od.away_spec{iC}.freqs, od.away_spec{iC}.sts_vals, 'blue');
%                     q0 = xline(away_sts_lg_peak, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     hold on;
%                     % Plot high gamma peak
%                     lf = find(od.away_spec{iC}.freqs >= hg(1), 1, 'first');
%                     rf = find(od.away_spec{iC}.freqs <= hg(2), 1, 'last');
%                     pks = findpeaks(od.away_spec{iC}.sts_vals(lf:rf));
%                     max_val = max(od.away_spec{iC}.sts_vals);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_spec{iC}.sts_vals(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                        [~,p1] = max(od.away_spec{iC}.sts_vals(lf:rf));
%                     end
%                     away_sts_hg_peak = od.away_spec{iC}.freqs(lf+p1-1);
%                     q0 = xline(away_sts_hg_peak, '--blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';      
%                     this_legend = {};
%                     flag_leg = true;
%                     this_legend{length(this_legend)+1} = sprintf ...
%                         ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
%                             away_sts_lg_peak, away_sts_hg_peak);
%                 else
% %                     close all;
% %                     continue;
%                 end
%                 if ~od.away_lfr_spec{iC}.flag_nansts && od.away_lfr_spec{iC}.spk_count >= 100
%                     % Plot low gamma peak
%                     lf = find(od.away_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
%                     rf = find(od.away_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
%                     pks = findpeaks(od.away_lfr_spec{iC}.sts_vals(lf:rf));
%                     max_val = max(od.away_lfr_spec{iC}.sts_vals);
%                      if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_lfr_spec{iC}.sts_vals(lf:rf));
%                     end
%                     away_lfr_sts_lg_peak = od.away_lfr_spec{iC}.freqs(lf+p1-1);
%                     plot(od.away_lfr_spec{iC}.freqs, od.away_lfr_spec{iC}.sts_vals, 'Color', 'red');
%                     q0 = xline(away_lfr_sts_lg_peak, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     % Plot high gamma peak
%                     lf = find(od.away_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
%                     rf = find(od.away_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
%                     pks = findpeaks(od.away_lfr_spec{iC}.sts_vals(lf:rf));
%                     max_val = max(od.away_lfr_spec{iC}.sts_vals);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_lfr_spec{iC}.sts_vals(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_lfr_spec{iC}.sts_vals(lf:rf));
%                     end
%                     away_lfr_sts_hg_peak = od.away_lfr_spec{iC}.freqs(lf+p1-1);
%                     q0 = xline(away_lfr_sts_hg_peak, '--red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                         ('LFR trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
%                             away_lfr_sts_lg_peak, away_lfr_sts_hg_peak);   
%                 else
% %                     close all;
% %                     continue;
%                 end
%                 if ~od.away_hfr_spec{iC}.flag_nansts && od.away_hfr_spec{iC}.spk_count >= 100
%                     % Plot low gamma peak
%                     lf = find(od.away_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
%                     rf = find(od.away_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
%                     pks = findpeaks(od.away_hfr_spec{iC}.sts_vals(lf:rf));
%                     max_val = max(od.away_hfr_spec{iC}.sts_vals);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_hfr_spec{iC}.sts_vals(lf:rf));
%                     end
%                     away_hfr_sts_lg_peak = od.away_hfr_spec{iC}.freqs(lf+p1-1);
%                     plot(od.away_hfr_spec{iC}.freqs, od.away_hfr_spec{iC}.sts_vals, 'Color', 'green');
%                     q0 = xline(away_hfr_sts_lg_peak, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     % Plot high gamma peak
%                     lf = find(od.away_hfr_spec{iC}.freqs >= hg(1), 1, 'first');
%                     rf = find(od.away_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
%                     pks = findpeaks(od.away_hfr_spec{iC}.sts_vals(lf:rf));
%                     max_val = max(od.away_hfr_spec{iC}.sts_vals);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_hfr_spec{iC}.sts_vals(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                        [~,p1] = max(od.away_hfr_spec{iC}.sts_vals(lf:rf));
%                     end
%                     away_hfr_sts_hg_peak = od.away_hfr_spec{iC}.freqs(lf+p1-1);
%                     q0 = xline(away_hfr_sts_hg_peak, '--green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                         ('HFR trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
%                             away_hfr_sts_lg_peak, away_hfr_sts_hg_peak); 
%                 else
% %                     close all;
% %                     continue;
%                 end
%                 if flag_leg
%                     legend(this_legend, 'Location', 'northwest')
%                     xlabel('Freqs')
%                     title('Away Reward STS');
%                 end
%                 
%                 flag_leg = false;
%                 subplot(3,3,9)
%                 if ~od.away_spec{iC}.flag_nanppc && od.away_spec{iC}.spk_count >= 100
%                     % Plot low gamma peak
%                     lf = find(od.away_spec{iC}.freqs >= lg(1), 1, 'first');
%                     rf = find(od.away_spec{iC}.freqs <= lg(2), 1, 'last');
%                     pks = findpeaks(od.away_spec{iC}.ppc(lf:rf));
%                     max_val = max(od.away_spec{iC}.ppc);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_spec{iC}.ppc(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_spec{iC}.ppc(lf:rf));
%                     end
%                     away_ppc_lg_peak = od.away_spec{iC}.freqs(lf+p1-1);
%                     plot(od.away_spec{iC}.freqs, od.away_spec{iC}.ppc, 'blue');
%                     q0 = xline(away_ppc_lg_peak, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     hold on;
%                     % Plot high gamma peak
%                     lf = find(od.away_spec{iC}.freqs >= hg(1), 1, 'first');
%                     rf = find(od.away_spec{iC}.freqs <= hg(2), 1, 'last');
%                     pks = findpeaks(od.away_spec{iC}.ppc(lf:rf));
%                     max_val = max(od.away_spec{iC}.ppc);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_spec{iC}.ppc(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                        [~,p1] = max(od.away_spec{iC}.ppc(lf:rf));
%                     end
%                     away_ppc_hg_peak = od.away_spec{iC}.freqs(lf+p1-1);
%                     q0 = xline(away_ppc_hg_peak, 'blue');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend = {};
%                     flag_leg = true;
%                     this_legend{length(this_legend)+1} = sprintf ...
%                         ('All Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
%                             away_ppc_lg_peak, away_ppc_hg_peak); 
%                 else
% %                     close all;
% %                     continue;
%                 end
%                 if ~od.away_lfr_spec{iC}.flag_nanppc && od.away_lfr_spec{iC}.spk_count >= 100
%                     % Plot low gamma peak
%                     lf = find(od.away_lfr_spec{iC}.freqs >= lg(1), 1, 'first');
%                     rf = find(od.away_lfr_spec{iC}.freqs <= lg(2), 1, 'last');
%                     pks = findpeaks(od.away_lfr_spec{iC}.ppc(lf:rf));
%                     max_val = max(od.away_lfr_spec{iC}.ppc);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_lfr_spec{iC}.ppc(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_lfr_spec{iC}.ppc(lf:rf));
%                     end
%                     away_lfr_ppc_lg_peak = od.away_lfr_spec{iC}.freqs(lf+p1-1);
%                     plot(od.away_lfr_spec{iC}.freqs, od.away_lfr_spec{iC}.ppc, 'Color', 'red');
%                     q0 = xline(away_lfr_ppc_lg_peak, 'red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     % Plot high gamma peak
%                     lf = find(od.away_lfr_spec{iC}.freqs >= hg(1), 1, 'first');
%                     rf = find(od.away_lfr_spec{iC}.freqs <= hg(2), 1, 'last');
%                     pks = findpeaks(od.away_lfr_spec{iC}.ppc(lf:rf));
%                     max_val = max(od.away_lfr_spec{iC}.ppc);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_lfr_spec{iC}.ppc(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_lfr_spec{iC}.ppc(lf:rf));
%                     end
%                     away_lfr_ppc_hg_peak = od.away_lfr_spec{iC}.freqs(lf+p1-1);
%                     q0 = xline(away_lfr_ppc_hg_peak, '--red');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                         ('LFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
%                             away_lfr_ppc_lg_peak, away_lfr_ppc_hg_peak); 
%                 else
% %                     close all;
% %                     continue;
%                 end
%                 if ~od.away_hfr_spec{iC}.flag_nanppc && od.away_hfr_spec{iC}.spk_count >= 100
%                     % Plot low gamma peak
%                     lf = find(od.away_hfr_spec{iC}.freqs >= lg(1), 1, 'first');
%                     rf = find(od.away_hfr_spec{iC}.freqs <= lg(2), 1, 'last');
%                     pks = findpeaks(od.away_hfr_spec{iC}.ppc(lf:rf));
%                     max_val = max(od.away_hfr_spec{iC}.ppc);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_hfr_spec{iC}.ppc(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_hfr_spec{iC}.ppc(lf:rf));
%                     end
%                     away_hfr_ppc_lg_peak = od.away_hfr_spec{iC}.freqs(lf+p1-1);
%                     plot(od.away_hfr_spec{iC}.freqs, od.away_hfr_spec{iC}.ppc, 'Color', 'green');
%                     q0 = xline(away_hfr_ppc_lg_peak, 'green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     % Plot high gamma peak
%                     lf = find(od.away_hfr_spec{iC}.freqs >=hg(1), 1, 'first');
%                     rf = find(od.away_hfr_spec{iC}.freqs <= hg(2), 1, 'last');
%                     pks = findpeaks(od.away_hfr_spec{iC}.ppc(lf:rf));
%                     max_val = max(od.away_hfr_spec{iC}.ppc);
%                     if ~isempty(pks.loc) 
%                         [pk1, loc1] = max(od.away_hfr_spec{iC}.ppc(lf+pks.loc-1));
%                         if pk1 > max_val*pk_thresh
%                             p1 = pks.loc(loc1);
%                         else
% %                            close all;
% %                            continue;
%                         end
%                     else
% %                         close all;
% %                         continue;
%                         [~,p1] = max(od.away_hfr_spec{iC}.ppc(lf:rf));
%                     end
%                     away_hfr_ppc_hg_peak = od.away_hfr_spec{iC}.freqs(lf+p1-1);
%                     q0 = xline(away_hfr_ppc_hg_peak, '--green');
%                     q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
%                     this_legend{length(this_legend)+1} = sprintf ...
%                         ('HFR Trials lg peak: %2.2f Hz, hg peak: %2.2f Hz ',...
%                             away_hfr_ppc_lg_peak, away_hfr_ppc_hg_peak);
%                 else
% %                     close all;
% %                     continue;
%                 end
%                 if flag_leg
%                     legend(this_legend, 'Location', 'northeast')
%                     xlabel('Freqs')
%                     title('Away Reward PPC');
%                 end
%             else
% %                 close all;
% %                 continue;
%             end
            
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
            if od.cell_type(iC) == 1
                o_name = cat(2, o_prefix, '_MSN');
            else
                 o_name = cat(2, o_prefix, '_FSI');
            end
            
           WriteFig(fig,o_name,1);
           close all;  
        end
    end
end
