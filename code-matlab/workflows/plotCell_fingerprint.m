cd('D:\RandomVstrAnalysis\combined_results_final');
rats = {'R117','R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
    
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        lf = find(od.S1.freqs >= 40, 1, 'first');
        rf = find(od.S1.freqs <= 80, 1, 'last');
        
        
        msn_labels = od.S1.label(od.S1.cell_type == 1);
        for iC  = 1:length(od.S1.msn_res)
            o_prefix = extractBefore(msn_labels(iC),'.t');
            % Skip cells with less than 200 spikes in either of the set of
            % trials
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end

            fig = figure('WindowState', 'maximized');
            
            subplot(5,1,1);
            [~,p1] = max(od.S1.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_l10');
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
            
            subplot(5,1,2);
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
            
            subplot(5,1,3);
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
            
            subplot(5,1,4);
            plot(od.S1.msn_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S1.msn_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S1.msn_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(5,1,5);
            lfr_sc = num2str(od.S1.msn_res(iC).scount_ptile(1));
            hfr_sc = num2str(od.S1.msn_res(iC).scount_ptile(2));
            lfr_fr = cat(2,num2str(od.S1.msn_res(iC).ptile_mfrs(1,1)), ...
                ' Hz - ',num2str(od.S1.msn_res(iC).ptile_mfrs(1,2)), ' Hz');
            hfr_fr = cat(2,num2str(od.S1.msn_res(iC).ptile_mfrs(2,1)), ...
                ' Hz - ',num2str(od.S1.msn_res(iC).ptile_mfrs(2,2)), ' Hz');
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
            
            subplot(5,1,1);
            [~,p1] = max(od.S1.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S1.freqs(lf+p1-1);
            hfr_max = od.S1.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_l10');
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
            
            subplot(5,1,2);
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
            
            subplot(5,1,3);
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
            
            subplot(5,1,4);
            plot(od.S1.fsi_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S1.fsi_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S1.fsi_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(5,1,5);
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
            
            subplot(5,1,1);
            [~,p1] = max(od.S2.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_l10');
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
            
            subplot(5,1,2);
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
            
            subplot(5,1,3);
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
            
            subplot(5,1,4);
            plot(od.S2.msn_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S2.msn_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S2.msn_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(5,1,5);
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
            
            subplot(5,1,1);
            [~,p1] = max(od.S2.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).mtsts_ptile(2,lf:rf));
            lfr_max = od.S2.freqs(lf+p1-1);
            hfr_max = od.S2.freqs(lf+p2-1);
            % Indicate category of cell in suffix
            if hfr_max - lfr_max > 20
                o_prefix = cat(2, o_prefix{1}, '_g20');
            elseif hfr_max - lfr_max < -20
                o_prefix = cat(2, o_prefix{1}, '_gn20');
            elseif hfr_max - lfr_max > 10
                o_prefix = cat(2, o_prefix{1}, '_g10');
            elseif hfr_max - lfr_max < -10
                o_prefix = cat(2, o_prefix{1}, '_gn10');
            else
                o_prefix = cat(2, o_prefix{1}, '_l10');
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
            
            subplot(5,1,2);
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
            
            subplot(5,1,3);
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
            
            subplot(5,1,4);
            plot(od.S2.fsi_res(iC).sta_ptile(1,:), 'blue')
            hold on;
            plot(od.S2.fsi_res(iC).sta_ptile(2,:), 'red')
            xlim([0, length(od.S2.fsi_res(iC).sta_ptile(1,:))+1])
            set(gca, 'XTickLabel', [-500:100:500])
            title('STA');
            
            subplot(5,1,5);
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
