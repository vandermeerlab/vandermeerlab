% Generate distribution of difference of frequency of peak SFc(in the range
% 40-80 hZ) for HFR and LFR trials
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/combined_results/');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*full.mat');
    ofiles = dir(searchString);
    
    msn_near_sfc_difs = [];
    msn_away_sfc_difs = [];
    fsi_near_usfc_difs = [];
    fsi_near_sfc_difs = [];
    fsi_away_usfc_difs = [];
    fsi_away_sfc_difs = [];
    
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        fq1 = find(od.S1.freqs >= 40, 1, 'first');
        fq2 = find(od.S1.freqs <= 80, 1, 'last');
        for iC  = 1:length(od.S1.msn_res)
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.msn_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S1.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.msn_res(iC).mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            msn_near_sfc_difs = [msn_near_sfc_difs, od.S1.freqs(fq1+id_hf) - od.S1.freqs(fq1+id_lf)];
        end
        
        for iC  = 1:length(od.S1.fsi_res)
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S1.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).mtsts_ptile(2,:);
            usfc_lf = od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).unsampled_mtsts_ptile(1,:);
            usfc_hf = od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).unsampled_mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            [~,uid_lf] = max(usfc_lf(fq1:fq2));
            [~,uid_hf] = max(usfc_hf(fq1:fq2));
            fsi_near_sfc_difs = [fsi_near_sfc_difs, od.S1.freqs(fq1+id_hf) - od.S1.freqs(fq1+id_lf)];
            fsi_near_usfc_difs = [fsi_near_usfc_difs, od.S1.freqs(fq1+uid_hf) - od.S1.freqs(fq1+uid_lf)];

        end
        
        fq1 = find(od.S2.freqs >= 40, 1, 'first');
        fq2 = find(od.S2.freqs <= 80, 1, 'last');
        for iC  = 1:length(od.S2.msn_res)
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S2.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S2.msn_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S2.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.msn_res(iC).mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            msn_away_sfc_difs = [msn_away_sfc_difs, od.S2.freqs(fq1+id_hf) - od.S2.freqs(fq1+id_lf)];
        end
        
        for iC  = 1:length(od.S2.fsi_res)
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S2.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).mtsts_ptile(2,:);
            usfc_lf = od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(1,:) ./ od.S2.fsi_res(iC).unsampled_mtsts_ptile(1,:);
            usfc_hf = od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).unsampled_mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf(fq1:fq2));
            [~,id_hf] = max(sfc_hf(fq1:fq2));
            [~,uid_lf] = max(usfc_lf(fq1:fq2));
            [~,uid_hf] = max(usfc_hf(fq1:fq2));
            fsi_away_sfc_difs = [fsi_away_sfc_difs, od.S2.freqs(fq1+id_hf) - od.S2.freqs(fq1+id_lf)];
            fsi_away_usfc_difs = [fsi_away_usfc_difs, od.S2.freqs(fq1+uid_hf) - od.S2.freqs(fq1+uid_lf)];

        end
    end
end
%%
figure
subplot(2,3,1)
histogram(msn_near_sfc_difs,-47.5:5:47.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_sfc_difs,-47.5:5:47.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Near Trials SFC (MSN: Red, FSI: Sampled, Green');

subplot(2,3,2)
histogram(msn_near_sfc_difs,-47.5:5:47.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_usfc_difs,-47.5:5:47.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Near Trials SFC (MSN: Red, FSI: Unsampled, Green');

subplot(2,3,3)
histogram(fsi_near_sfc_difs,-47.5:5:47.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_usfc_difs,-47.5:5:47.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Near Trials SFC (Sampled FSI : Red, Unsampled FSI : Green');

subplot(2,3,4)
histogram(msn_away_sfc_difs,-47.5:5:47.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_sfc_difs,-47.5:5:47.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Away Trials SFC (MSN: Red, FSI: Sampled, Green');

subplot(2,3,5)
histogram(msn_away_sfc_difs,-47.5:5:47.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_usfc_difs,-47.5:5:47.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Away Trials SFC (MSN: Red, FSI: Unsampled, Green');

subplot(2,3,6)
histogram(fsi_away_sfc_difs,-47.5:5:47.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_usfc_difs,-47.5:5:47.5, 'FaceColor','green','FaceAlpha',0.4)
hold on
title('Away Trials SFC (Sampled FSI : Red, Unsampled FSI : Green');


