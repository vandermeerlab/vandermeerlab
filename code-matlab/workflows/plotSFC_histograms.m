% Plot distribution of peak SFC frequency (in the range lfq Hz - hfq Hz) for
% all cases: fsi-sampled, fsi-unsampled, msn
lfq = 0;
hfq = 100;
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/Results/combined_results_optimal/');
rats = {'R117','R119','R131','R132'};
    msn_near_sfc_lf = [];
    msn_near_sfc_hf = [];
    msn_away_sfc_lf = [];
    msn_away_sfc_hf = [];
    fsi_near_usfc_lf = [];
    fsi_near_usfc_hf = [];
    fsi_near_sfc_lf = [];
    fsi_near_sfc_hf = [];
    fsi_away_usfc_lf = [];
    fsi_away_usfc_hf = [];
    fsi_away_sfc_lf = [];
    fsi_away_sfc_hf = [];
    
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
    
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        fq1 = find(od.S1.freqs >= lfq, 1, 'first');
        fq2 = find(od.S1.freqs <= hfq, 1, 'last');
        for iC  = 1:length(od.S1.msn_res)
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.msn_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S1.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.msn_res(iC).mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            msn_near_sfc_lf = [msn_near_sfc_lf, od.S1.freqs(fq1+id_lf-1)];
            msn_near_sfc_hf = [msn_near_sfc_hf, od.S1.freqs(fq1+id_hf-1)];
        end
        
        for iC  = 1:length(od.S1.fsi_res)
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S1.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).mtsts_ptile(2,:);
            usfc_lf = od.S1.fsi_res(iC).us_sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).us_mtsts_ptile(1,:);
            usfc_hf = od.S1.fsi_res(iC).us_sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).us_mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            [~,uid_lf] = max(usfc_lf(fq1:fq2));
            [~,uid_hf] = max(usfc_hf(fq1:fq2));
            fsi_near_sfc_lf = [fsi_near_sfc_lf, od.S1.freqs(fq1+id_lf-1)];
            fsi_near_sfc_hf = [fsi_near_sfc_hf, od.S1.freqs(fq1+id_hf-1)];
            fsi_near_usfc_lf = [fsi_near_usfc_lf, od.S1.freqs(fq1+uid_lf-1)];
            fsi_near_usfc_hf = [fsi_near_usfc_hf, od.S1.freqs(fq1+uid_hf-1)];
        end
        
        fq1 = find(od.S2.freqs >= lfq, 1, 'first');
        fq2 = find(od.S2.freqs <= hfq, 1, 'last');
        for iC  = 1:length(od.S2.msn_res)
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S2.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S2.msn_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S2.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.msn_res(iC).mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            msn_away_sfc_lf = [msn_away_sfc_lf, od.S2.freqs(fq1+id_lf-1)];
            msn_away_sfc_hf = [msn_away_sfc_hf, od.S2.freqs(fq1+id_hf-1)];
        end
        
        for iC  = 1:length(od.S2.fsi_res)
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S2.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).mtsts_ptile(2,:);
            usfc_lf = od.S2.fsi_res(iC).us_sta_mtspec_ptile(1,:) ./ od.S2.fsi_res(iC).us_mtsts_ptile(1,:);
            usfc_hf = od.S2.fsi_res(iC).us_sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).us_mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            [~,uid_lf] = max(usfc_lf(fq1:fq2));
            [~,uid_hf] = max(usfc_hf(fq1:fq2));
            fsi_away_sfc_lf = [fsi_away_sfc_lf, od.S2.freqs(fq1+id_lf-1)];
            fsi_away_sfc_hf = [fsi_away_sfc_hf, od.S2.freqs(fq1+id_hf-1)];
            fsi_away_usfc_lf = [fsi_away_usfc_lf, od.S2.freqs(fq1+uid_lf-1)];
            fsi_away_usfc_hf = [fsi_away_usfc_hf, od.S2.freqs(fq1+uid_hf-1)];
        end
    end
end
%%
close all;
figure
subplot(2,3,1)
histogram(msn_near_sfc_lf,-4.5:5:100.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(msn_near_sfc_hf,-4.5:5:100.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Near Trials MSN SFC (LFR : Green, HFR : Red)');

subplot(2,3,2)
histogram(fsi_near_sfc_lf,-4.5:5:100.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_sfc_hf,-4.5:5:100.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Near Trials FSI SFC (1000 times Sampled)');

subplot(2,3,3)
histogram(fsi_near_usfc_lf,-4.5:5:100.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_usfc_hf,-4.5:5:100.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Near Trials FSI SFC (Unsampled)');

subplot(2,3,4)
histogram(msn_away_sfc_lf,-4.5:5:100.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(msn_away_sfc_hf,-4.5:5:100.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Away Trials MSN SFC (LFR : Green, HFR : Red)');

subplot(2,3,5)
histogram(fsi_away_sfc_lf,-4.5:5:100.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_sfc_hf,-4.5:5:100.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Away Trials FSI SFC (1000 times Sampled)');

subplot(2,3,6)
histogram(fsi_away_usfc_lf,-4.5:5:100.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_usfc_hf,-4.5:5:100.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Away Trials FSI SFC (Unsampled)');

