% Generate distribution of peak SFC frequency for all cases
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/combined_results/');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*full.mat');
    ofiles = dir(searchString);
    
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
    
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
%         lf = find(od.S1.freqs >= 40, 1, 'first');
%         rf = find(od.S1.freqs <= 80, 1, 'last');
        for iC  = 1:length(od.S1.msn_res)
            if od.S1.msn_res(iC).scount_ptile(1) < 200 || od.S1.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.msn_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S1.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.msn_res(iC).mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf);
            [~,id_hf] = max(sfc_hf);
            msn_near_sfc_lf = [msn_near_sfc_lf, od.S1.freqs(id_lf)];
            msn_near_sfc_hf = [msn_near_sfc_hf, od.S1.freqs(id_hf)];
        end
        
        for iC  = 1:length(od.S1.fsi_res)
            if od.S1.fsi_res(iC).scount_ptile(1) < 200 || od.S1.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S1.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).mtsts_ptile(2,:);
            usfc_lf = od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).unsampled_mtsts_ptile(1,:);
            usfc_hf = od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile(2,:) ./ od.S1.fsi_res(iC).unsampled_mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf);
            [~,id_hf] = max(sfc_hf);
            [~,uid_lf] = max(usfc_lf);
            [~,uid_hf] = max(usfc_hf);
            fsi_near_sfc_lf = [fsi_near_sfc_lf, od.S1.freqs(id_lf)];
            fsi_near_sfc_hf = [fsi_near_sfc_hf, od.S1.freqs(id_hf)];
            fsi_near_usfc_lf = [fsi_near_usfc_lf, od.S1.freqs(uid_lf)];
            fsi_near_usfc_hf = [fsi_near_usfc_hf, od.S1.freqs(uid_hf)];
        end
        
        for iC  = 1:length(od.S2.msn_res)
            if od.S2.msn_res(iC).scount_ptile(1) < 200 || od.S2.msn_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S2.msn_res(iC).sta_mtspec_ptile(1,:) ./ od.S2.msn_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S2.msn_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.msn_res(iC).mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf);
            [~,id_hf] = max(sfc_hf);
            msn_away_sfc_lf = [msn_away_sfc_lf, od.S2.freqs(id_lf)];
            msn_away_sfc_hf = [msn_away_sfc_hf, od.S2.freqs(id_hf)];
        end
        
        for iC  = 1:length(od.S2.fsi_res)
            if od.S2.fsi_res(iC).scount_ptile(1) < 200 || od.S2.fsi_res(iC).scount_ptile(2) < 200
                continue
            end
            sfc_lf = od.S1.fsi_res(iC).sta_mtspec_ptile(1,:) ./ od.S1.fsi_res(iC).mtsts_ptile(1,:);
            sfc_hf = od.S2.fsi_res(iC).sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).mtsts_ptile(2,:);
            usfc_lf = od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(1,:) ./ od.S2.fsi_res(iC).unsampled_mtsts_ptile(1,:);
            usfc_hf = od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile(2,:) ./ od.S2.fsi_res(iC).unsampled_mtsts_ptile(2,:);
            [~,id_lf] = max(sfc_lf);
            [~,id_hf] = max(sfc_hf);
            [~,uid_lf] = max(usfc_lf);
            [~,uid_hf] = max(usfc_hf);
            fsi_away_sfc_lf = [fsi_away_sfc_lf, od.S2.freqs(id_lf)];
            fsi_away_sfc_hf = [fsi_away_sfc_hf, od.S2.freqs(id_hf)];
            fsi_away_usfc_lf = [fsi_away_usfc_lf, od.S2.freqs(uid_lf)];
            fsi_away_usfc_hf = [fsi_away_usfc_hf, od.S2.freqs(uid_hf)];
        end
    end
end
%%
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

