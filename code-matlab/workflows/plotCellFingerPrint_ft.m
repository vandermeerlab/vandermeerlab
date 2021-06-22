cd('D:\RandomVstrAnalysis\combined_results');
rats = {'R117','R119','R131','R132'};
lfq = 30;
hfq = 100;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    
    for jdx = 1%:length(ofiles)
        load(ofiles(jdx).name);
        dummy = 2;
        num_cells = length(od.label);
        for iC = 1:num_cells    
            % Plot On Track stuff first
            subplot(3,3,1)
            plot(od.onTrack_spec{iC}.sta_time, od.onTrack_spec{iC}.sta_vals);
            title(sprintf('On Track STA, %d spikes', od.onTrack_spec{iC}.spk_count));
            
            if ~od.onTrack_spec{iC}.flag_nansts
                subplot(3,3,2)
                plot(od.onTrack_spec{iC}.freqs, od.onTrack_spec{iC}.sts_vals);
                xlabel('Freqs')
                title('On Track STS');
            end
            
            if ~od.onTrack_spec{iC}.flag_nanppc
                subplot(3,3,3)
                plot(od.onTrack_spec{iC}.freqs, od.onTrack_spec{iC}.ppc);
                xlabel('Freqs')
                title('On Track PPC');
            end
            
            % Plot Near Reward stuff next
            if ~od.near_spec{iC}.flag_zeroSpikes
                
                subplot(3,3,4)
                plot(od.near_spec{iC}.sta_time, od.near_spec{iC}.sta_vals);
                hold on;
                this_legend = {};
                this_legend{1} = sprintf('Near All Spikes: %d',od.near_spec{iC}.spk_count);
                if ~od.near_lfr_spec{iC}.flag_zeroSpikes
                    plot(od.near_spec{iC}.sta_time, od.near_lfr_spec{iC}.sta_vals, 'Color', 'red');
                    this_legend{2} = sprintf('Near LFR spikes: %d',od.near_lfr_spec{iC}.spk_count);
                end
                if ~od.near_hfr_spec{iC}.flag_zeroSpikes
                    plot(od.near_spec{iC}.sta_time, od.near_hfr_spec{iC}.sta_vals, 'Color', 'green');
                    this_legend{3} = sprintf('Near HFR spikes: %d',od.near_hfr_spec{iC}.spk_count);
                end
                legend(this_legend)
                title('Near Reward STA');
                
                
                
                
            end

            
            if ~od.onTrack_spec{iC}.flag_nansts
                subplot(3,3,2)
                plot(od.onTrack_spec{iC}.freqs, od.onTrack_spec{iC}.sts_vals);
                xlabel('Freqs')
                title('On Track STS');
            end
            
            if ~od.onTrack_spec{iC}.flag_nanppc
                subplot(3,3,3)
                plot(od.onTrack_spec{iC}.freqs, od.onTrack_spec{iC}.ppc);
                xlabel('Freqs')
                title('On Track PPC');
            end
            
            
            
        end
        
        msn_labels = od.S1.label(od.S1.cell_type == 1);
        for iC  = 1:length(od.S1.msn_res)
            o_prefix = extractBefore(msn_labels(iC),'.t');
            % Skip cells with less than 200 spikes in either of the set of
            % trials
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end

            fig = figure('WindowState', 'maximized');
            
            subplot(6,1,1);
            [~,p1] = max(od.S1.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_sts_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_sts_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_sts_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_sts_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_sts_l10');
            end
                     
            plot(od.S1.freqs, 10*log(od.S1.msn_res(iC).mtsts_ptile(1,:)), 'blue');
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S1.freqs, 10*log(od.S1.msn_res(iC).mtsts_ptile(2,:)), 'red');
            xline(hfr_max, 'red');
            [m1,~] = max(od.S1.msn_res(iC).mtsts_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STS');
            
            subplot(6,1,2);
            [~,p1] = max(od.S1.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            plot(od.S1.freqs, 10*log(od.S1.msn_res(iC).sta_mtspec_ptile(1,:)), 'blue')
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S1.freqs, 10*log(od.S1.msn_res(iC).sta_mtspec_ptile(2,:)), 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(od.S1.msn_res(iC).sta_mtspec_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STA Spectrum');
            
            subplot(6,1,3);
            [~,p1] = max(od.S1.msn_res(iC).psd(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).psd(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            plot(od.S1.freqs, 10*log(od.S1.msn_res(iC).psd(1,:)), 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S1.freqs, 10*log(od.S1.msn_res(iC).psd(2,:)), 'red')
            xline(hfr_max, 'red')
            [m1,~] = max(od.S1.msn_res(iC).psd(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('LFP Spectrum');
            
            subplot(6,1,4);
            plot(od.S1.msn_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S1.msn_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S1.msn_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(6,1,5);
            msn_sfc_lf = od.S1.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.msn_res(iC).mtsts_ptile(1,:);
            msn_sfc_hf = od.S1.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.msn_res(iC).mtsts_ptile(2,:);
            [~,p1] = max(msn_sfc_lf(lf:rf));
            [~,p2] = max(msn_sfc_hf(lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix, '_sfc_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix, '_sfc_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix, '_sfc_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix, '_sfc_gn10');
            else
                o_prefix = cat(2, o_prefix, '_sfc_l10');
            end
            
            plot(od.S1.freqs, msn_sfc_lf, 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S1.freqs, msn_sfc_hf, 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(msn_sfc_lf);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
%             text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
%             text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
%             text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('SFC');
                 
            %Calculate mean_lfr and mean_hfr
            lfr_mask = od.S1.msn_res(iC).mfr <= od.S1.msn_res(iC).ptile_mfrs(1,2);
            hfr_mask = ~(lfr_mask);
            nzfr_mask = (od.S1.msn_res(iC).mfr ~= 0);
            lfr_mask = lfr_mask & nzfr_mask;
            hfr_mask = hfr_mask & nzfr_mask;
            trial_length = od.S1.trial_ends - od.S1.trial_starts;
            lfr_time = sum(trial_length(lfr_mask));
            hfr_time = sum(trial_length(hfr_mask));
            mean_lfr = od.S1.msn_res(iC).scount_ptile(1)/lfr_time;
            mean_hfr = od.S1.msn_res(iC).scount_ptile(2)/hfr_time;
            if mean_lfr > od.S1.msn_res(iC).ptile_mfrs(1,2) || mean_hfr <= od.S1.msn_res(iC).ptile_mfrs(1,2)
                o_prefix = cat(2, o_prefix, '_fr_ERROR');
            end
             % Indicate category of cell in suffix
            if mean_hfr - mean_lfr <= 3
                o_prefix = cat(2, o_prefix, '_fr_l3');
            elseif mean_hfr - mean_lfr <= 6
                o_prefix = cat(2, o_prefix, '_fr_l6');
            else
                o_prefix = cat(2, o_prefix, '_fr_g6');
            end
            
            subplot(6,1,6);
            lfr_sc = num2str(od.S1.msn_res(iC).scount_ptile(1));
            hfr_sc = num2str(od.S1.msn_res(iC).scount_ptile(2));
            lfr_fr = cat(2,num2str(od.S1.msn_res(iC).ptile_mfrs(1,1)), ...
                ' Hz - ',num2str(od.S1.msn_res(iC).ptile_mfrs(1,2)), ' Hz');
            hfr_fr = cat(2,num2str(od.S1.msn_res(iC).ptile_mfrs(2,1)), ...
                ' Hz - ',num2str(od.S1.msn_res(iC).ptile_mfrs(2,2)), ' Hz');
          
            text(0.05, 1.1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(0.35, 1.1, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(0.65, 1.1, dif_txt, 'Color', 'black', 'FontSize', 12);           
            text(0.2,0.8,cat(2,'Spike Count for low firing rate trials : ',...
                lfr_sc), 'FontSize',14);
            text(0.2,0.5,cat(2,'Spike Count for high firing rate trials : ',...
                hfr_sc), 'FontSize', 14);
            text(0.5,0.8,cat(2,'Firing rate range for low firing rate trials : ',...
                lfr_fr), 'FontSize', 14);
            text(0.5,0.5,cat(2,'Firing rate range for high firing rate trials : ',...
                hfr_fr), 'FontSize', 14);
            set(gca,'xtick', []);
            set(gca,'ytick', []);   
            dummy  = 22;
            o_name = cat(2, o_prefix, '_Near_MSN');
            WriteFig(fig,o_name,1);
            close all;
        end
        
        fsi_labels = od.S1.label(od.S1.cell_type == 2);
        for iC  = 1:length(od.S1.fsi_res)
            o_prefix = extractBefore(fsi_labels(iC),'.t');
            % Skip cells with less than 200 spikes in either of the set of
            % trials
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end

            fig = figure('WindowState', 'maximized');
            
            subplot(6,1,1);
            [~,p1] = max(od.S1.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_sts_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_sts_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_sts_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_sts_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_sts_l10');
            end
                     
            plot(od.S1.freqs, 10*log(od.S1.fsi_res(iC).mtsts_ptile(1,:)), 'blue');
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S1.freqs, 10*log(od.S1.fsi_res(iC).mtsts_ptile(2,:)), 'red');
            xline(hfr_max, 'red');
            [m1,~] = max(od.S1.fsi_res(iC).mtsts_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STS');
            
            subplot(6,1,2);
            [~,p1] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            plot(od.S1.freqs, 10*log(od.S1.fsi_res(iC).sta_mtspec_ptile(1,:)), 'blue')
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S1.freqs, 10*log(od.S1.fsi_res(iC).sta_mtspec_ptile(2,:)), 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STA Spectrum');
            
            subplot(6,1,3);
            [~,p1] = max(od.S1.fsi_res(iC).psd(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).psd(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            plot(od.S1.freqs, 10*log(od.S1.fsi_res(iC).psd(1,:)), 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S1.freqs, 10*log(od.S1.fsi_res(iC).psd(2,:)), 'red')
            xline(hfr_max, 'red')
            [m1,~] = max(od.S1.fsi_res(iC).psd(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('LFP Spectrum');
            
            subplot(6,1,4);
            plot(od.S1.fsi_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S1.fsi_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S1.fsi_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(6,1,5);
            fsi_sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            fsi_sfc_hf = od.S1.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).mtsts_ptile(2,:);
            [~,p1] = max(fsi_sfc_lf(lf:rf));
            [~,p2] = max(fsi_sfc_hf(lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix, '_sfc_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix, '_sfc_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix, '_sfc_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix, '_sfc_gn10');
            else
                o_prefix = cat(2, o_prefix, '_sfc_l10');
            end
            
            plot(od.S1.freqs, fsi_sfc_lf, 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S1.freqs, fsi_sfc_hf, 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(fsi_sfc_lf);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
%             text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
%             text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
%             text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('SFC');
            
            %Calculate mean_lfr and mean_hfr
            lfr_mask = od.S1.fsi_res(iC).mfr <= od.S1.fsi_res(iC).ptile_mfrs(1,2);
            hfr_mask = ~(lfr_mask);
            nzfr_mask = (od.S1.fsi_res(iC).mfr ~= 0);
            lfr_mask = lfr_mask & nzfr_mask;
            hfr_mask = hfr_mask & nzfr_mask;
            trial_length = od.S1.trial_ends - od.S1.trial_starts;
            lfr_time = sum(trial_length(lfr_mask));
            hfr_time = sum(trial_length(hfr_mask));
            mean_lfr = od.S1.fsi_res(iC).scount_ptile(1)/lfr_time;
            mean_hfr = od.S1.fsi_res(iC).scount_ptile(2)/hfr_time;
            if mean_lfr > od.S1.fsi_res(iC).ptile_mfrs(1,2) || mean_hfr <= od.S1.fsi_res(iC).ptile_mfrs(1,2)
                o_prefix = cat(2, o_prefix, '_fr_ERROR');
            end
             % Indicate category of cell in suffix
            if mean_hfr - mean_lfr <= 3
                o_prefix = cat(2, o_prefix, '_fr_l3');
            elseif mean_hfr - mean_lfr <= 6
                o_prefix = cat(2, o_prefix, '_fr_l6');
            else
                o_prefix = cat(2, o_prefix, '_fr_g6');
            end
            
            subplot(6,1,6);
            text(0.05, 1.1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(0.35, 1.1, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(0.65, 1.1, dif_txt, 'Color', 'black', 'FontSize', 12); 
            lfr_sc = num2str(od.S1.fsi_res(iC).scount_ptile(1));
            hfr_sc = num2str(od.S1.fsi_res(iC).scount_ptile(2));
            lfr_fr = cat(2,num2str(od.S1.fsi_res(iC).ptile_mfrs(1,1)), ...
                ' Hz - ',num2str(od.S1.fsi_res(iC).ptile_mfrs(1,2)), ' Hz');
            hfr_fr = cat(2,num2str(od.S1.fsi_res(iC).ptile_mfrs(2,1)), ...
                ' Hz - ',num2str(od.S1.fsi_res(iC).ptile_mfrs(2,2)), ' Hz');
            text(0.2,0.8,cat(2,'Spike Count for low firing rate trials : ',...
                lfr_sc), 'FontSize',14);
            text(0.2,0.5,cat(2,'Spike Count for high firing rate trials : ',...
                hfr_sc), 'FontSize', 14);
            text(0.5,0.8,cat(2,'Firing rate range for low firing rate trials : ',...
                lfr_fr), 'FontSize', 14);
            text(0.5,0.5,cat(2,'Firing rate range for high firing rate trials : ',...
                hfr_fr), 'FontSize', 14);
            set(gca,'xtick', []);
            set(gca,'ytick', []);   
            dummy  = 22;
            o_name = cat(2, o_prefix, '_Near_FSI');
            WriteFig(fig,o_name,1);
            close all;
        end
    
        msn_labels = od.S2.label(od.S2.cell_type == 1);
        for iC  = 1:length(od.S2.msn_res)
            o_prefix = extractBefore(msn_labels(iC),'.t');
            % Skip cells with less than 200 spikes in either of the set of
            % trials
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end

            fig = figure('WindowState', 'maximized');
            
            subplot(6,1,1);
            [~,p1] = max(od.S2.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_sts_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_sts_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_sts_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_sts_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_sts_l10');
            end
                     
            plot(od.S2.freqs, 10*log(od.S2.msn_res(iC).mtsts_ptile(1,:)), 'blue');
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S2.freqs, 10*log(od.S2.msn_res(iC).mtsts_ptile(2,:)), 'red');
            xline(hfr_max, 'red');
            [m1,~] = max(od.S2.msn_res(iC).mtsts_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STS');
            
            subplot(6,1,2);
            [~,p1] = max(od.S2.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            plot(od.S2.freqs, 10*log(od.S2.msn_res(iC).sta_mtspec_ptile(1,:)), 'blue')
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S2.freqs, 10*log(od.S2.msn_res(iC).sta_mtspec_ptile(2,:)), 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(od.S2.msn_res(iC).sta_mtspec_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STA Spectrum');
            
            subplot(6,1,3);
            [~,p1] = max(od.S2.msn_res(iC).psd(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).psd(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            plot(od.S2.freqs, 10*log(od.S2.msn_res(iC).psd(1,:)), 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S2.freqs, 10*log(od.S2.msn_res(iC).psd(2,:)), 'red')
            xline(hfr_max, 'red')
            [m1,~] = max(od.S2.msn_res(iC).psd(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('LFP Spectrum');
            
            subplot(6,1,4);
            plot(od.S2.msn_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S2.msn_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S2.msn_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(6,1,5);
            msn_sfc_lf = od.S2.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S2.msn_res(iC).mtsts_ptile(1,:);
            msn_sfc_hf = od.S2.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.msn_res(iC).mtsts_ptile(2,:);
            [~,p1] = max(msn_sfc_lf(lf:rf));
            [~,p2] = max(msn_sfc_hf(lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix, '_sfc_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix, '_sfc_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix, '_sfc_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix, '_sfc_gn10');
            else
                o_prefix = cat(2, o_prefix, '_sfc_l10');
            end
            
            plot(od.S2.freqs, msn_sfc_lf, 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S2.freqs, msn_sfc_hf, 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(msn_sfc_lf);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
%             text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
%             text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
%             text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('SFC');
            
            %Calculate mean_lfr and mean_hfr
            lfr_mask = od.S2.msn_res(iC).mfr <= od.S2.msn_res(iC).ptile_mfrs(1,2);
            hfr_mask = ~(lfr_mask);
            nzfr_mask = (od.S2.msn_res(iC).mfr ~= 0);
            lfr_mask = lfr_mask & nzfr_mask;
            hfr_mask = hfr_mask & nzfr_mask;
            trial_length = od.S2.trial_ends - od.S2.trial_starts;
            lfr_time = sum(trial_length(lfr_mask));
            hfr_time = sum(trial_length(hfr_mask));
            mean_lfr = od.S2.msn_res(iC).scount_ptile(1)/lfr_time;
            mean_hfr = od.S2.msn_res(iC).scount_ptile(2)/hfr_time;
            if mean_lfr > od.S2.msn_res(iC).ptile_mfrs(1,2) || mean_hfr <= od.S2.msn_res(iC).ptile_mfrs(1,2)
                o_prefix = cat(2, o_prefix, '_fr_ERROR');
            end
             % Indicate category of cell in suffix
            if mean_hfr - mean_lfr <= 3
                o_prefix = cat(2, o_prefix, '_fr_l3');
            elseif mean_hfr - mean_lfr <= 6
                o_prefix = cat(2, o_prefix, '_fr_l6');
            else
                o_prefix = cat(2, o_prefix, '_fr_g6');
            end
            
            subplot(6,1,6);
            text(0.05, 1.1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(0.35, 1.1, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(0.65, 1.1, dif_txt, 'Color', 'black', 'FontSize', 12); 
            lfr_sc = num2str(od.S2.msn_res(iC).scount_ptile(1));
            hfr_sc = num2str(od.S2.msn_res(iC).scount_ptile(2));
            lfr_fr = cat(2,num2str(od.S2.msn_res(iC).ptile_mfrs(1,1)), ...
                ' Hz - ',num2str(od.S2.msn_res(iC).ptile_mfrs(1,2)), ' Hz');
            hfr_fr = cat(2,num2str(od.S2.msn_res(iC).ptile_mfrs(2,1)), ...
                ' Hz - ',num2str(od.S2.msn_res(iC).ptile_mfrs(2,2)), ' Hz');
            text(0.2,0.8,cat(2,'Spike Count for low firing rate trials : ',...
                lfr_sc), 'FontSize',14);
            text(0.2,0.5,cat(2,'Spike Count for high firing rate trials : ',...
                hfr_sc), 'FontSize', 14);
            text(0.5,0.8,cat(2,'Firing rate range for low firing rate trials : ',...
                lfr_fr), 'FontSize', 14);
            text(0.5,0.5,cat(2,'Firing rate range for high firing rate trials : ',...
                hfr_fr), 'FontSize', 14);
            set(gca,'xtick', []);
            set(gca,'ytick', []);   
            dummy  = 22;
            o_name = cat(2, o_prefix, '_Away_MSN');
            WriteFig(fig,o_name,1);
            close all;
        end
        
        fsi_labels = od.S2.label(od.S2.cell_type == 2);
        for iC  = 1:length(od.S2.fsi_res)
            o_prefix = extractBefore(fsi_labels(iC),'.t');
            % Skip cells with less than 200 spikes in either of the set of
            % trials
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end

            fig = figure('WindowState', 'maximized');
            
            subplot(6,1,1);
            [~,p1] = max(od.S2.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_sts_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_sts_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_sts_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_sts_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_sts_l10');
            end
                     
            plot(od.S2.freqs, 10*log(od.S2.fsi_res(iC).mtsts_ptile(1,:)), 'blue');
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S2.freqs, 10*log(od.S2.fsi_res(iC).mtsts_ptile(2,:)), 'red');
            xline(hfr_max, 'red');
            [m1,~] = max(od.S2.fsi_res(iC).mtsts_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STS');
            
            subplot(6,1,2);
            [~,p1] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            plot(od.S2.freqs, 10*log(od.S2.fsi_res(iC).sta_mtspec_ptile(1,:)), 'blue')
            xline(lfr_max, 'blue');
            hold on;
            plot(od.S2.freqs, 10*log(od.S2.fsi_res(iC).sta_mtspec_ptile(2,:)), 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('STA Spectrum');
            
            subplot(6,1,3);
            [~,p1] = max(od.S2.fsi_res(iC).psd(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).psd(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            plot(od.S2.freqs, 10*log(od.S2.fsi_res(iC).psd(1,:)), 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S2.freqs, 10*log(od.S2.fsi_res(iC).psd(2,:)), 'red')
            xline(hfr_max, 'red')
            [m1,~] = max(od.S2.fsi_res(iC).psd(1,:));
            m1 = 10*log(m1);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
            text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('LFP Spectrum');
            
            subplot(6,1,4);
            plot(od.S2.fsi_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S2.fsi_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S2.fsi_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(6,1,5);
            fsi_sfc_lf = od.S2.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S2.fsi_res(iC).mtsts_ptile(1,:);
            fsi_sfc_hf = od.S2.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).mtsts_ptile(2,:);
            [~,p1] = max(fsi_sfc_lf(lf:rf));
            [~,p2] = max(fsi_sfc_hf(lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix, '_sfc_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix, '_sfc_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix, '_sfc_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix, '_sfc_gn10');
            else
                o_prefix = cat(2, o_prefix, '_sfc_l10');
            end
            
            plot(od.S2.freqs, fsi_sfc_lf, 'blue')
            xline(lfr_max, 'blue')
            hold on;
            plot(od.S2.freqs, fsi_sfc_hf, 'red')
            xline(hfr_max, 'red');
            [m1,~] = max(fsi_sfc_lf);
            lfr_txt = cat(2, 'Peak Freq for Low Firing Rate Trials : ', ...
                num2str(lfr_max), ' Hz');
            hfr_txt = cat(2, 'Peak Freq for High Firing Rate Trials : ', ...
                num2str(hfr_max), ' Hz');
            dif_txt = cat(2, 'HFR Peak - LFR Peak : ', ...
                num2str(hfr_max - lfr_max), ' Hz');
%             text(80, m1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
%             text(80, m1 - 10, hfr_txt, 'Color', 'red', 'FontSize', 12);
%             text(80, m1 - 20, dif_txt, 'Color', 'black', 'FontSize', 12);
            title('SFC');
            
            %Calculate mean_lfr and mean_hfr
            lfr_mask = od.S2.fsi_res(iC).mfr <= od.S2.fsi_res(iC).ptile_mfrs(1,2);
            hfr_mask = ~(lfr_mask);
            nzfr_mask = (od.S2.fsi_res(iC).mfr ~= 0);
            lfr_mask = lfr_mask & nzfr_mask;
            hfr_mask = hfr_mask & nzfr_mask;
            trial_length = od.S2.trial_ends - od.S2.trial_starts;
            lfr_time = sum(trial_length(lfr_mask));
            hfr_time = sum(trial_length(hfr_mask));
            mean_lfr = od.S2.fsi_res(iC).scount_ptile(1)/lfr_time;
            mean_hfr = od.S2.fsi_res(iC).scount_ptile(2)/hfr_time;
            if mean_lfr > od.S2.fsi_res(iC).ptile_mfrs(1,2) || mean_hfr <= od.S2.fsi_res(iC).ptile_mfrs(1,2)
                o_prefix = cat(2, o_prefix, '_fr_ERROR');
            end
             % Indicate category of cell in suffix
            if mean_hfr - mean_lfr <= 3
                o_prefix = cat(2, o_prefix, '_fr_l3');
            elseif mean_hfr - mean_lfr <= 6
                o_prefix = cat(2, o_prefix, '_fr_l6');
            else
                o_prefix = cat(2, o_prefix, '_fr_g6');
            end
            
            subplot(6,1,6);
            text(0.05, 1.1, lfr_txt, 'Color', 'blue', 'FontSize', 12);
            text(0.35, 1.1, hfr_txt, 'Color', 'red', 'FontSize', 12);
            text(0.65, 1.1, dif_txt, 'Color', 'black', 'FontSize', 12); 
            lfr_sc = num2str(od.S2.fsi_res(iC).scount_ptile(1));
            hfr_sc = num2str(od.S2.fsi_res(iC).scount_ptile(2));
            lfr_fr = cat(2,num2str(od.S2.fsi_res(iC).ptile_mfrs(1,1)), ...
                ' Hz - ',num2str(od.S2.fsi_res(iC).ptile_mfrs(1,2)), ' Hz');
            hfr_fr = cat(2,num2str(od.S2.fsi_res(iC).ptile_mfrs(2,1)), ...
                ' Hz - ',num2str(od.S2.fsi_res(iC).ptile_mfrs(2,2)), ' Hz');
            text(0.2,0.8,cat(2,'Spike Count for low firing rate trials : ',...
                lfr_sc), 'FontSize',14);
            text(0.2,0.5,cat(2,'Spike Count for high firing rate trials : ',...
                hfr_sc), 'FontSize', 14);
            text(0.5,0.8,cat(2,'Firing rate range for low firing rate trials : ',...
                lfr_fr), 'FontSize', 14);
            text(0.5,0.5,cat(2,'Firing rate range for high firing rate trials : ',...
                hfr_fr), 'FontSize', 14);
            set(gca,'xtick', []);
            set(gca,'ytick', []);   
            dummy  = 22;
            o_name = cat(2, o_prefix, '_Away_FSI');
            WriteFig(fig,o_name,1);
            close all;
        end
    end
end
