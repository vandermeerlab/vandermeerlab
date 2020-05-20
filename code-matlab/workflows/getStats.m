cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/sts_v4');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
    msn_near_sta_difs = [];
    msn_away_sta_difs = [];
    msn_near_sts_difs = [];
    msn_away_sts_difs = [];
    fsi_near_sta_difs = [];
    fsi_away_sta_difs = [];
    fsi_near_sts_difs = [];
    fsi_away_sts_difs = [];
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        lf = find(od.S1.freqs >= 40, 1, 'first');
        rf = find(od.S1.freqs <= 80, 1, 'last');
        for iC  = 1:length(od.S1.cell_type)
            [~,p1] = max(od.S1.spec_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S1.spec_res(iC).mtsts_ptile(2,lf:rf));
            if od.S1.cell_type(iC) == 1
                msn_near_sts_difs = [msn_near_sts_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            else
                fsi_near_sts_difs = [fsi_near_sts_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            end
            [~,p1] = max(od.S1.spec_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p2] = max(od.S1.spec_res(iC).sta_mtspec_ptile(2,lf:rf));
            if od.S1.cell_type(iC) == 1
                msn_near_sta_difs = [msn_near_sta_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            else
                fsi_near_sta_difs = [fsi_near_sta_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            end
        end
        lf = find(od.S2.freqs >= 40, 1, 'first');
        rf = find(od.S2.freqs <= 80, 1, 'last');
        for iC  = 1:length(od.S2.cell_type)
            [~,p1] = max(od.S2.spec_res(iC).mtsts_ptile(1,lf:rf));
            [~,p2] = max(od.S2.spec_res(iC).mtsts_ptile(2,lf:rf));
            if od.S2.cell_type(iC) == 1
                msn_away_sts_difs = [msn_away_sts_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            else
                fsi_away_sts_difs = [fsi_away_sts_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            end
            [~,p1] = max(od.S2.spec_res(iC).sta_mtspec_ptile(1,lf:rf));
            [~,p2] = max(od.S2.spec_res(iC).sta_mtspec_ptile(2,lf:rf));
            if od.S2.cell_type(iC) == 1
                msn_away_sta_difs = [msn_away_sta_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            else
                fsi_away_sta_difs = [fsi_away_sta_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            end
        end
    end
end
subplot(2,2,1)
histogram(msn_away_sts_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_sts_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Away trials STS');
subplot(2,2,3)
histogram(msn_near_sts_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_sts_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Near trials STS');
subplot(2,2,2)
histogram(msn_away_sta_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_sta_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Away trials STA spec');
subplot(2,2,4)
histogram(msn_near_sta_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_sta_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Near trials STA spec');