% Generate scatter plots of (MFR of High FR trials - MFR of Low FR Trials) vs various quantities

lfq = 0;
hfq = 100;
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/Results/combined_results_optimal/');
rats = {'R117','R119','R131','R132'};

msn_near_fr_difs = [];
msn_away_fr_difs = [];
msn_near_sta_difs = [];
msn_away_sta_difs = [];
msn_near_sts_difs = [];
msn_away_sts_difs = [];
msn_near_lfr_sfc = [];
msn_near_hfr_sfc = [];
msn_away_lfr_sfc = [];
msn_away_hfr_sfc = [];
msn_near_psd_difs = [];
msn_away_psd_difs = [];


fsi_near_fr_difs = [];
fsi_away_fr_difs = [];
fsi_near_sta_difs = [];
fsi_away_sta_difs = [];
fsi_near_sts_difs = [];
fsi_away_sts_difs = [];
fsi_near_lfr_sfc = [];
fsi_near_hfr_sfc = [];
fsi_away_lfr_sfc = [];
fsi_away_hfr_sfc = [];
fsi_near_psd_difs = [];
fsi_away_psd_difs = [];

for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);   
    for jdx = 1:length(ofiles)

        load(ofiles(jdx).name);
        lf = find(od.S1.freqs >= lfq, 1, 'first');
        rf = find(od.S1.freqs <= hfq, 1, 'last');
       
        % Skip cells that have less than 200 spikes in either of the two
        % trial sets
        for iC  = 1:length(od.S1.msn_res)
            od.S1.msn_res(iC).sfc = od.S1.msn_res(iC).sta_mtspec_ptile ./ od.S1.msn_res(iC).mtsts_ptile;
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end
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
                continue
            end
            msn_near_fr_difs = [msn_near_fr_difs mean_hfr-mean_lfr];
          
            [~,p1] = max(od.S1.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            [lfr_sfc,~] = max(od.S1.msn_res(iC).sfc(1,lf:rf));
            [hfr_sfc,~] = max(od.S1.msn_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S1.msn_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S1.msn_res(iC).psd(2,lf:rf));
            
            msn_near_sts_difs = [msn_near_sts_difs (od.S1.freqs(lf+p2-1)-od.S1.freqs(lf+p1-1))];
            msn_near_sta_difs = [msn_near_sta_difs (od.S1.freqs(lf+p4-1)-od.S1.freqs(lf+p3-1))];     
            msn_near_lfr_sfc = [msn_near_lfr_sfc lfr_sfc];
            msn_near_hfr_sfc = [msn_near_hfr_sfc hfr_sfc];      
            msn_near_psd_difs = [msn_near_psd_difs (od.S1.freqs(lf+p8-1)-od.S1.freqs(lf+p7-1))];
        end
        
         for iC  = 1:length(od.S1.fsi_res)
            od.S1.fsi_res(iC).sfc = od.S1.fsi_res(iC).sta_mtspec_ptile ./ od.S1.fsi_res(iC).mtsts_ptile;
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            
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
                continue
            end
            fsi_near_fr_difs = [fsi_near_fr_difs mean_hfr-mean_lfr];
            
            [~,p1] = max(od.S1.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            [lfr_sfc,~] = max(od.S1.fsi_res(iC).sfc(1,lf:rf));
            [hfr_sfc,~] = max(od.S1.fsi_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S1.fsi_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S1.fsi_res(iC).psd(2,lf:rf));

            fsi_near_sts_difs = [fsi_near_sts_difs (od.S1.freqs(lf+p2-1)-od.S1.freqs(lf+p1-1))];
            fsi_near_sta_difs = [fsi_near_sta_difs (od.S1.freqs(lf+p4-1)-od.S1.freqs(lf+p3-1))];
            fsi_near_lfr_sfc = [fsi_near_lfr_sfc lfr_sfc];
            fsi_near_hfr_sfc = [fsi_near_hfr_sfc hfr_sfc]; 
            fsi_near_psd_difs = [fsi_near_psd_difs (od.S1.freqs(lf+p8-1)-od.S1.freqs(lf+p7-1))];
        end
        
        lf = find(od.S2.freqs >= lfq, 1, 'first');
        rf = find(od.S2.freqs <= hfq, 1, 'last');
        
        for iC  = 1:length(od.S2.msn_res)
            od.S2.msn_res(iC).sfc = od.S2.msn_res(iC).sta_mtspec_ptile ./ od.S2.msn_res(iC).mtsts_ptile;
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            
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
                continue
            end
            msn_away_fr_difs = [msn_away_fr_difs mean_hfr-mean_lfr];
            
            [~,p1] = max(od.S2.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            [lfr_sfc,~] = max(od.S2.msn_res(iC).sfc(1,lf:rf));
            [hfr_sfc,~] = max(od.S2.msn_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S2.msn_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S2.msn_res(iC).psd(2,lf:rf));
            
            msn_away_sts_difs = [msn_away_sts_difs (od.S2.freqs(lf+p2-1)-od.S2.freqs(lf+p1-1))];
            msn_away_sta_difs = [msn_away_sta_difs (od.S2.freqs(lf+p4-1)-od.S2.freqs(lf+p3-1))];
            msn_away_lfr_sfc = [msn_away_lfr_sfc lfr_sfc];
            msn_away_hfr_sfc = [msn_away_hfr_sfc hfr_sfc];
            msn_away_psd_difs = [msn_away_psd_difs (od.S2.freqs(lf+p8-1)-od.S2.freqs(lf+p7-1))];
        end
        
        for iC  = 1:length(od.S2.fsi_res)
            od.S2.fsi_res(iC).sfc = od.S2.fsi_res(iC).sta_mtspec_ptile ./ od.S2.fsi_res(iC).mtsts_ptile;
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            
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
                continue
            end
            fsi_away_fr_difs = [fsi_away_fr_difs mean_hfr-mean_lfr];
            
            [~,p1] = max(od.S2.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            [lfr_sfc,~] = max(od.S2.fsi_res(iC).sfc(1,lf:rf));
            [hfr_sfc,~] = max(od.S2.fsi_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S2.fsi_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S2.fsi_res(iC).psd(2,lf:rf));

            fsi_away_sts_difs = [fsi_away_sts_difs (od.S2.freqs(lf+p2-1)-od.S2.freqs(lf+p1-1))];
            fsi_away_sta_difs = [fsi_away_sta_difs (od.S2.freqs(lf+p4-1)-od.S2.freqs(lf+p3-1))];
            fsi_away_lfr_sfc = [fsi_away_lfr_sfc lfr_sfc];
            fsi_away_hfr_sfc = [fsi_away_hfr_sfc hfr_sfc];
            fsi_away_psd_difs = [fsi_away_psd_difs (od.S2.freqs(lf+p8-1)-od.S2.freqs(lf+p7-1))];
        end
        

        
    end
end
msn_near_sfc = msn_near_lfr_sfc./2 + msn_near_hfr_sfc./2;
fsi_near_sfc = fsi_near_lfr_sfc./2 + fsi_near_hfr_sfc./2;
msn_away_sfc = msn_away_lfr_sfc./2 + msn_away_hfr_sfc./2;
fsi_away_sfc = fsi_away_lfr_sfc./2 + fsi_away_hfr_sfc./2;


close all;

figure;
subplot(2,4,1);
scatter(msn_near_fr_difs, msn_near_sts_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near MFR Dif')
ylabel('MSN Near STS Dif')
subplot(2,4,2);
scatter(msn_near_fr_difs, msn_near_sta_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near MFR Dif')
ylabel('MSN Near STA Dif')
subplot(2,4,3);
scatter(msn_near_fr_difs, msn_near_psd_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near MFR Dif')
ylabel('MSN Near PSD Dif')
subplot(2,4,4);
scatter(msn_near_fr_difs, msn_near_sfc, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near MFR Dif')
ylabel('MSN Near mean SFC')


subplot(2,4,5);
scatter(fsi_near_fr_difs, fsi_near_sts_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near MFR Dif')
ylabel('FSI Near STS Dif')
subplot(2,4,6);
scatter(fsi_near_fr_difs, fsi_near_sta_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near MFR Dif')
ylabel('FSI Near STA Dif')
subplot(2,4,7);
scatter(fsi_near_fr_difs, fsi_near_psd_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near MFR Dif')
ylabel('FSI Near PSD Dif')
subplot(2,4,8);
scatter(fsi_near_fr_difs, fsi_near_sfc, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near MFR Dif')
ylabel('FSI Near mean SFC')

figure;
subplot(2,4,1);
scatter(msn_away_fr_difs, msn_away_sts_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away MFR Dif')
ylabel('MSN Away STS Dif')
subplot(2,4,2);
scatter(msn_away_fr_difs, msn_away_sta_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away MFR Dif')
ylabel('MSN Away STA Dif')
subplot(2,4,3);
scatter(msn_away_fr_difs, msn_away_psd_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away MFR Dif')
ylabel('MSN Away PSD Dif')
subplot(2,4,4);
scatter(msn_away_fr_difs, msn_away_sfc, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away MFR Dif')
ylabel('MSN Away mean SFC')


subplot(2,4,5);
scatter(fsi_away_fr_difs, fsi_away_sts_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away MFR Dif')
ylabel('FSI Away STS Dif')
subplot(2,4,6);
scatter(fsi_away_fr_difs, fsi_away_sta_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away MFR Dif')
ylabel('FSI Away STA Dif')
subplot(2,4,7);
scatter(fsi_away_fr_difs, fsi_away_psd_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away MFR Dif')
ylabel('FSI Away PSD Dif')
subplot(2,4,8);
scatter(fsi_away_fr_difs, fsi_away_sfc, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away MFR Dif')
ylabel('FSI Away mean SFC')

figure
subplot(2,1,1)
for i = 1:length(msn_away_fr_difs)
plot([msn_away_fr_difs(i),msn_away_fr_difs(i)],[msn_away_lfr_sfc(i),msn_away_hfr_sfc(i)],'Color', 'r');
hold on
end
xlabel('MSN Away FR Dif')
ylabel('MSN Away SFC (low-high)')
subplot(2,1,2)
for i = 1:length(fsi_away_fr_difs)
plot([fsi_away_fr_difs(i),fsi_away_fr_difs(i)],[fsi_away_lfr_sfc(i),fsi_away_hfr_sfc(i)],'Color', 'g');
hold on
end
xlabel('FSI Away FR Dif')
ylabel('FSI Away SFC (low-high)')

figure
subplot(2,1,1)
for i = 1:length(msn_near_fr_difs)
plot([msn_near_fr_difs(i),msn_near_fr_difs(i)],[msn_near_lfr_sfc(i),msn_near_hfr_sfc(i)],'Color', 'r');
hold on
end
xlabel('MSN Near FR Dif')
ylabel('MSN Near SFC (low-high)')
subplot(2,1,2)
for i = 1:length(fsi_near_fr_difs)
plot([fsi_near_fr_difs(i),fsi_near_fr_difs(i)],[fsi_near_lfr_sfc(i),fsi_near_hfr_sfc(i)],'Color', 'g');
hold on
end
xlabel('FSI Near FR Dif')
ylabel('FSI Near SFC (low-high)')

