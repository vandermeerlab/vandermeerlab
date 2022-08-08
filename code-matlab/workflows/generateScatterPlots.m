cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/Results/combined_results_optimal/');
rats = {'R117','R119','R131','R132'};

msn_near_sta_difs = [];
msn_away_sta_difs = [];
msn_near_sts_difs = [];
msn_away_sts_difs = [];
msn_near_sfc_difs = [];
msn_away_sfc_difs = [];
msn_near_psd_difs = [];
msn_away_psd_difs = [];

fsi_near_sta_difs = [];
fsi_away_sta_difs = [];
fsi_near_sts_difs = [];
fsi_away_sts_difs = [];
fsi_near_sfc_difs = [];
fsi_away_sfc_difs = [];
fsi_near_psd_difs = [];
fsi_away_psd_difs = [];


for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
%     
    for jdx = 1:length(ofiles)

        load(ofiles(jdx).name);
        lf = find(od.S1.freqs >= 40, 1, 'first');
        rf = find(od.S1.freqs <= 80, 1, 'last');
       
        
        % Skip cells that have less than 200 spikes in either of the two
        % trial sets
        for iC  = 1:length(od.S1.msn_res)
            od.S1.msn_res(iC).sfc = od.S1.msn_res(iC).sta_mtspec_ptile ./ od.S1.msn_res(iC).mtsts_ptile;
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S1.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S1.msn_res(iC).sfc(1,lf:rf));
            [~,p6] = max(od.S1.msn_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S1.msn_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S1.msn_res(iC).psd(2,lf:rf));
            
            msn_near_sts_difs = [msn_near_sts_difs (od.S1.freqs(lf+p2-1)-od.S1.freqs(lf+p1-1))];
            msn_near_sta_difs = [msn_near_sta_difs (od.S1.freqs(lf+p4-1)-od.S1.freqs(lf+p3-1))];
            msn_near_sfc_difs = [msn_near_sfc_difs (od.S1.freqs(lf+p6-1)-od.S1.freqs(lf+p5-1))];
            msn_near_psd_difs = [msn_near_psd_difs (od.S1.freqs(lf+p8-1)-od.S1.freqs(lf+p7-1))];
        end
        
         for iC  = 1:length(od.S1.fsi_res)
            od.S1.fsi_res(iC).sfc = od.S1.fsi_res(iC).sta_mtspec_ptile ./ od.S1.fsi_res(iC).mtsts_ptile;
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S1.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S1.fsi_res(iC).sfc(1,lf:rf));
            [~,p6] = max(od.S1.fsi_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S1.fsi_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S1.fsi_res(iC).psd(2,lf:rf));

            fsi_near_sts_difs = [fsi_near_sts_difs (od.S1.freqs(lf+p2-1)-od.S1.freqs(lf+p1-1))];
            fsi_near_sta_difs = [fsi_near_sta_difs (od.S1.freqs(lf+p4-1)-od.S1.freqs(lf+p3-1))];
            fsi_near_sfc_difs = [fsi_near_sfc_difs (od.S1.freqs(lf+p6-1)-od.S1.freqs(lf+p5-1))];
            fsi_near_psd_difs = [fsi_near_psd_difs (od.S1.freqs(lf+p8-1)-od.S1.freqs(lf+p7-1))];
        end
        
        lf = find(od.S2.freqs >= 40, 1, 'first');
        rf = find(od.S2.freqs <= 80, 1, 'last');
        
        for iC  = 1:length(od.S2.msn_res)
            od.S2.msn_res(iC).sfc = od.S2.msn_res(iC).sta_mtspec_ptile ./ od.S2.msn_res(iC).mtsts_ptile;
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S2.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S2.msn_res(iC).sfc(1,lf:rf));
            [~,p6] = max(od.S2.msn_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S2.msn_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S2.msn_res(iC).psd(2,lf:rf));
            
            msn_away_sts_difs = [msn_away_sts_difs (od.S2.freqs(lf+p2-1)-od.S2.freqs(lf+p1-1))];
            msn_away_sta_difs = [msn_away_sta_difs (od.S2.freqs(lf+p4-1)-od.S2.freqs(lf+p3-1))];
            msn_away_sfc_difs = [msn_away_sfc_difs (od.S2.freqs(lf+p6-1)-od.S2.freqs(lf+p5-1))];
            msn_away_psd_difs = [msn_away_psd_difs (od.S2.freqs(lf+p8-1)-od.S2.freqs(lf+p7-1))];
        end
        
        for iC  = 1:length(od.S2.fsi_res)
            od.S2.fsi_res(iC).sfc = od.S2.fsi_res(iC).sta_mtspec_ptile ./ od.S2.fsi_res(iC).mtsts_ptile;
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S2.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S2.fsi_res(iC).sfc(1,lf:rf));
            [~,p6] = max(od.S2.fsi_res(iC).sfc(2,lf:rf));
            [~,p7] = max(od.S2.fsi_res(iC).psd(1,lf:rf));
            [~,p8] = max(od.S2.fsi_res(iC).psd(2,lf:rf));

            fsi_away_sts_difs = [fsi_away_sts_difs (od.S2.freqs(lf+p2-1)-od.S2.freqs(lf+p1-1))];
            fsi_away_sta_difs = [fsi_away_sta_difs (od.S2.freqs(lf+p4-1)-od.S2.freqs(lf+p3-1))];
            fsi_away_sfc_difs = [fsi_away_sfc_difs (od.S2.freqs(lf+p6-1)-od.S2.freqs(lf+p5-1))];
            fsi_away_psd_difs = [fsi_away_psd_difs (od.S2.freqs(lf+p8-1)-od.S2.freqs(lf+p7-1))];
        end
        

        
    end
end
close all;
figure;
subplot(2,3,1);
scatter(msn_near_sts_difs, msn_near_sta_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near STS Dif')
ylabel('MSN Near STA Dif')
subplot(2,3,2);
scatter(msn_near_sts_difs, msn_near_sfc_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near STS Dif')
ylabel('MSN Near SFC Dif')
subplot(2,3,3);
scatter(msn_near_sts_difs, msn_near_psd_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near STS Dif')
ylabel('MSN Near PSD Dif')
subplot(2,3,4);
scatter(msn_near_sta_difs, msn_near_sfc_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near STA Dif')
ylabel('MSN Near SFC Dif')
subplot(2,3,5);
scatter(msn_near_sta_difs, msn_near_psd_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near STA Dif')
ylabel('MSN Near PSD Dif')
subplot(2,3,6);
scatter(msn_near_psd_difs, msn_near_sfc_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Near PSD Dif')
ylabel('MSN Near SFC Dif')

figure;
subplot(2,3,1);
scatter(fsi_near_sts_difs, fsi_near_sta_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near STS Dif')
ylabel('FSI Near STA Dif')
subplot(2,3,2);
scatter(fsi_near_sts_difs, fsi_near_sfc_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near STS Dif')
ylabel('FSI Near SFC Dif')
subplot(2,3,3);
scatter(fsi_near_sts_difs, fsi_near_psd_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near STS Dif')
ylabel('FSI Near PSD Dif')
subplot(2,3,4);
scatter(fsi_near_sta_difs, fsi_near_sfc_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near STA Dif')
ylabel('FSI Near SFC Dif')
subplot(2,3,5);
scatter(fsi_near_sta_difs, fsi_near_psd_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near STA Dif')
ylabel('FSI Near PSD Dif')
subplot(2,3,6);
scatter(fsi_near_psd_difs, fsi_near_sfc_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Near PSD Dif')
ylabel('FSI Near SFC Dif')

figure;
subplot(2,3,1);
scatter(msn_away_sts_difs, msn_away_sta_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away STS Dif')
ylabel('MSN Away STA Dif')
subplot(2,3,2);
scatter(msn_away_sts_difs, msn_away_sfc_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away STS Dif')
ylabel('MSN Away SFC Dif')
subplot(2,3,3);
scatter(msn_away_sts_difs, msn_away_psd_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away STS Dif')
ylabel('MSN Away PSD Dif')
subplot(2,3,4);
scatter(msn_away_sta_difs, msn_away_sfc_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away STA Dif')
ylabel('MSN Away SFC Dif')
subplot(2,3,5);
scatter(msn_away_sta_difs, msn_away_psd_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away STA Dif')
ylabel('MSN Away PSD Dif')
subplot(2,3,6);
scatter(msn_away_psd_difs, msn_away_sfc_difs, 'filled', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.4);
xlabel('MSN Away PSD Dif')
ylabel('MSN Away SFC Dif')

figure;
subplot(2,3,1);
scatter(fsi_away_sts_difs, fsi_away_sta_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away STS Dif')
ylabel('FSI Away STA Dif')
subplot(2,3,2);
scatter(fsi_away_sts_difs, fsi_away_sfc_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away STS Dif')
ylabel('FSI Away SFC Dif')
subplot(2,3,3);
scatter(fsi_away_sts_difs, fsi_away_psd_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away STS Dif')
ylabel('FSI Away PSD Dif')
subplot(2,3,4);
scatter(fsi_away_sta_difs, fsi_away_sfc_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away STA Dif')
ylabel('FSI Away SFC Dif')
subplot(2,3,5);
scatter(fsi_away_sta_difs, fsi_away_psd_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away STA Dif')
ylabel('FSI Away PSD Dif')
subplot(2,3,6);
scatter(fsi_away_psd_difs, fsi_away_sfc_difs, 'filled', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.4);
xlabel('FSI Away PSD Dif')
ylabel('FSI Away SFC Dif')
