% Generate distribution of differences between frequency
% (in the range lfq Hz - hfq Hz) with peak STA and peak STS values in low 
% firing vs high firing trials
lfq = 0;
hfq = 100;
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/Results/combined_results_optimal/');
rats = {'R117','R119','R131','R132'};
msn_near_sta_difs = [];
msn_away_sta_difs = [];
msn_near_sts_difs = [];
msn_away_sts_difs = [];
fsi_near_sta_difs = [];
fsi_away_sta_difs = [];
fsi_near_sts_difs = [];
fsi_away_sts_difs = [];
msn_near_psd_difs = [];
msn_away_psd_difs = [];
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
        for iC  = 1:length(od.S1.msn_res)
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S1.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.msn_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S1.msn_res(iC).psd(1,lf:rf));
            [~,p6] = max(od.S1.msn_res(iC).psd(2,lf:rf));
            msn_near_sts_difs = [msn_near_sts_difs (od.S1.freqs(lf+p2-1)-od.S1.freqs(lf+p1-1))];
            msn_near_sta_difs = [msn_near_sta_difs (od.S1.freqs(lf+p4-1)-od.S1.freqs(lf+p3-1))];
            msn_near_psd_difs = [msn_near_psd_difs (od.S1.freqs(lf+p5-1)-od.S1.freqs(lf+p6-1))];
        end
        for iC  = 1:length(od.S1.fsi_res)
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S1.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
%             [~,p5] = max(od.S1.fsi_res(iC).unsampled_mtsts_ptile(1,lf:rf));
%             [~,p6] = max(od.S1.fsi_res(iC).unsampled_mtsts_ptile(2,lf:rf));
%             [~,p7] = max(od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(1,lf:rf));
%             [~,p8] = max(od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(2,lf:rf));
            [~,p9] = max(od.S1.fsi_res(iC).psd(1,lf:rf));
            [~,p10] = max(od.S1.fsi_res(iC).psd(2,lf:rf));
            fsi_near_sts_difs = [fsi_near_sts_difs (od.S1.freqs(lf+p2-1)-od.S1.freqs(lf+p1-1))];
            fsi_near_sta_difs = [fsi_near_sta_difs (od.S1.freqs(lf+p4-1)-od.S1.freqs(lf+p3-1))];
%             fsi_near_usts_difs = [fsi_near_usts_difs (od.S1.freqs(lf+p6-1)-od.S1.freqs(lf+p5-1))];
%             fsi_near_usta_difs = [fsi_near_usta_difs (od.S1.freqs(lf+p8-1)-od.S1.freqs(lf+p7-1))];
            fsi_near_psd_difs = [fsi_near_psd_difs (od.S1.freqs(lf+p10-1)-od.S1.freqs(lf+p9-1))];
        end
        lf = find(od.S2.freqs >= lfq, 1, 'first');
        rf = find(od.S2.freqs <= hfq, 1, 'last');
        for iC  = 1:length(od.S2.msn_res)
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S2.msn_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.msn_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.msn_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.msn_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S2.msn_res(iC).psd(1,lf:rf));
            [~,p6] = max(od.S2.msn_res(iC).psd(2,lf:rf));
            msn_away_sts_difs = [msn_away_sts_difs (od.S2.freqs(lf+p2-1)-od.S2.freqs(lf+p1-1))];
            msn_away_sta_difs = [msn_away_sta_difs (od.S2.freqs(lf+p4-1)-od.S2.freqs(lf+p3-1))];
            msn_away_psd_difs = [msn_away_psd_difs (od.S2.freqs(lf+p5-1)-od.S2.freqs(lf+p6-1))];
        end
        for iC  = 1:length(od.S2.fsi_res)
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S2.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
%             [~,p5] = max(od.S2.fsi_res(iC).unsampled_mtsts_ptile(1,lf:rf));
%             [~,p6] = max(od.S2.fsi_res(iC).unsampled_mtsts_ptile(2,lf:rf));
%             [~,p7] = max(od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(1,lf:rf));
%             [~,p8] = max(od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(2,lf:rf));
            [~,p9] = max(od.S2.fsi_res(iC).psd(1,lf:rf));
            [~,p10] = max(od.S2.fsi_res(iC).psd(2,lf:rf));
            fsi_away_sts_difs = [fsi_away_sts_difs (od.S2.freqs(lf+p2-1)-od.S2.freqs(lf+p1-1))];
            fsi_away_sta_difs = [fsi_away_sta_difs (od.S2.freqs(lf+p4-1)-od.S2.freqs(lf+p3-1))];
%             fsi_away_usts_difs = [fsi_away_usts_difs (od.S2.freqs(lf+p6-1)-od.S2.freqs(lf+p5-1))];
%             fsi_away_usta_difs = [fsi_away_usta_difs (od.S2.freqs(lf+p8-1)-od.S2.freqs(lf+p7-1))];
            fsi_away_psd_difs = [fsi_away_psd_difs (od.S2.freqs(lf+p10-1)-od.S2.freqs(lf+p9-1))];
        end
    end
end
close all;
figure;
subplot(2,3,1)
histogram(msn_away_sts_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_sts_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Away trials STS');
subplot(2,3,4)
histogram(msn_near_sts_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_sts_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Near trials STS');
subplot(2,3,2)
histogram(msn_away_sta_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_sta_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Away trials STA spec');
subplot(2,3,5)
histogram(msn_near_sta_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_sta_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Near trials STA spec');
subplot(2,3,3)
histogram(msn_away_psd_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_psd_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Away trials average LFP PSD');
subplot(2,3,6)
histogram(msn_near_psd_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_psd_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Near trials average LFP PSD');