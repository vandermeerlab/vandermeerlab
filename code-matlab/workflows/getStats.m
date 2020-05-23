cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/sts_v4');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
    msn_near_lfr_scount = [];
    msn_away_lfr_scount = [];
    fsi_near_scount = [];
    fsi_away_scount = [];
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        for iC  = 1:length(od.S1.cell_type)
            if od.S1.cell_type(iC) == 1
                msn_near_scount = [msn_near_scount sum(od.S1.spec_res(iC).scount_ptile)];
            else
                fsi_near_scount = [fsi_near_scount sum(od.S1.spec_res(iC).scount_ptile)];
            end
        end
        for iC  = 1:length(od.S2.cell_type)
            if od.S2.cell_type(iC) == 1
                msn_away_scount = [msn_away_scount sum(od.S2.spec_res(iC).scount_ptile)];
            else
                fsi_away_scount = [fsi_away_scount sum(od.S2.spec_res(iC).scount_ptile)];
            end
        end
    end
end
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