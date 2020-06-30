cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/combined_results/');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*full.mat');
    ofiles = dir(searchString);
   
    msn_near_sta_difs = [];
    msn_away_sta_difs = [];
    msn_near_sts_difs = [];
    msn_away_sts_difs = [];
    fsi_near_sta_difs = [];
    fsi_away_sta_difs = [];
    fsi_near_sts_difs = [];
    fsi_away_sts_difs = [];

    fsi_near_usta_difs = [];
    fsi_away_usta_difs = [];
    fsi_near_usts_difs = [];
    fsi_away_usts_difs = [];
    
    msn_near_psd_difs = [];
    msn_away_psd_difs = [];
    fsi_near_psd_difs = [];
    fsi_away_psd_difs = [];
    
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        lf = find(od.S1.freqs >= 40, 1, 'first');
        rf = find(od.S1.freqs <= 80, 1, 'last');
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
            msn_near_sts_difs = [msn_near_sts_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            msn_near_sta_difs = [msn_near_sta_difs (od.S1.freqs(lf+p4)-od.S1.freqs(lf+p3))];
            msn_near_psd_difs = [msn_near_psd_difs (od.S1.freqs(lf+p5)-od.S1.freqs(lf+p6))];
        end
        for iC  = 1:length(od.S1.fsi_res)
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S1.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S1.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S1.fsi_res(iC).unsampled_mtsts_ptile(1,lf:rf));
            [~,p6] = max(od.S1.fsi_res(iC).unsampled_mtsts_ptile(2,lf:rf));
            [~,p7] = max(od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(1,lf:rf));
            [~,p8] = max(od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(2,lf:rf));
            [~,p9] = max(od.S1.fsi_res(iC).psd(1,lf:rf));
            [~,p10] = max(od.S1.fsi_res(iC).psd(2,lf:rf));
            fsi_near_sts_difs = [fsi_near_sts_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            fsi_near_sta_difs = [fsi_near_sta_difs (od.S1.freqs(lf+p4)-od.S1.freqs(lf+p3))];
            fsi_near_usts_difs = [fsi_near_usts_difs (od.S1.freqs(lf+p6)-od.S1.freqs(lf+p5))];
            fsi_near_usta_difs = [fsi_near_usta_difs (od.S1.freqs(lf+p8)-od.S1.freqs(lf+p7))];
            fsi_near_psd_difs = [fsi_near_psd_difs (od.S1.freqs(lf+p10)-od.S1.freqs(lf+p9))];
        end
        lf = find(od.S2.freqs >= 40, 1, 'first');
        rf = find(od.S2.freqs <= 80, 1, 'last');
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
            msn_away_sts_difs = [msn_away_sts_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            msn_away_sta_difs = [msn_away_sta_difs (od.S2.freqs(lf+p4)-od.S2.freqs(lf+p3))];
            msn_away_psd_difs = [msn_away_psd_difs (od.S2.freqs(lf+p5)-od.S2.freqs(lf+p6))];
        end
        for iC  = 1:length(od.S2.fsi_res)
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S2.fsi_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.fsi_res(iC).mtsts_ptile(2,lf:rf));
            [~,p3] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p4] = max(od.S2.fsi_res(iC).sta_mtspec_ptile(2,lf:rf));
            [~,p5] = max(od.S2.fsi_res(iC).unsampled_mtsts_ptile(1,lf:rf));
            [~,p6] = max(od.S2.fsi_res(iC).unsampled_mtsts_ptile(2,lf:rf));
            [~,p7] = max(od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(1,lf:rf));
            [~,p8] = max(od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(2,lf:rf));
            [~,p9] = max(od.S2.fsi_res(iC).psd(1,lf:rf));
            [~,p10] = max(od.S2.fsi_res(iC).psd(2,lf:rf));
            fsi_away_sts_difs = [fsi_away_sts_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            fsi_away_sta_difs = [fsi_away_sta_difs (od.S2.freqs(lf+p4)-od.S2.freqs(lf+p3))];
            fsi_away_usts_difs = [fsi_away_usts_difs (od.S2.freqs(lf+p6)-od.S2.freqs(lf+p5))];
            fsi_away_usta_difs = [fsi_away_usta_difs (od.S2.freqs(lf+p8)-od.S2.freqs(lf+p7))];
            fsi_away_psd_difs = [fsi_away_psd_difs (od.S2.freqs(lf+p10)-od.S2.freqs(lf+p9))];
        end
    end
end
% figure;
% subplot(2,2,1)
% histogram(msn_away_sts_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
% hold on
% histogram(fsi_away_sts_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
% title('Away trials STS');
% subplot(2,2,3)
% histogram(msn_near_sts_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
% hold on
% histogram(fsi_near_sts_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
% title('Near trials STS');
% subplot(2,2,2)
% histogram(msn_away_sta_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
% hold on
% histogram(fsi_away_sta_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
% title('Away trials STA spec');
% subplot(2,2,4)
% histogram(msn_near_sta_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
% hold on
% histogram(fsi_near_sta_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
% title('Near trials STA spec');