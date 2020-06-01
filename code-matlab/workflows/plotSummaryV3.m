cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/psd_v1');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*psd.mat');
    ofiles = dir(searchString);
    msn_near_psd_difs = [];
    msn_away_psd_difs = [];
    fsi_near_psd_difs = [];
    fsi_away_psd_difs = [];
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        lf = find(od.S1.freqs >= 40, 1, 'first');
        rf = find(od.S1.freqs <= 80, 1, 'last');
        
        % Hack to exclude sessions in which none of the MSNs have at least
        % 200 spikes in both the percentiles
        msn1 = find(od.S1.cell_type == 1);
        throw_msn1 = 0;
        for iC = 1:length(msn1)
            if od.S1.spec_res(msn1(iC)).scount_ptile(1) < 200 || od.S1.spec_res(msn1(iC)).scount_ptile(2) < 200
                throw_msn1 = throw_msn1 + 1;
            end
        end
        msn2 = find(od.S2.cell_type == 1);
        throw_msn2 = 0;
        for iC = 1:length(msn2)
            if od.S2.spec_res(msn2(iC)).scount_ptile(1) < 200 || od.S2.spec_res(msn2(iC)).scount_ptile(2) < 200
                throw_msn2 = throw_msn2 + 1;
            end
        end
        
        for iC  = 1:length(od.S1.cell_type)
            if throw_msn1 == length(msn1) || od.S1.spec_res(iC).scount_ptile(1) < 200 || od.S1.spec_res(iC).scount_ptile(2) < 200 
                continue
            end
            [~,p1] = max(od.S1.spec_res(iC).psd(1,lf:rf));
            [~,p2] = max(od.S1.spec_res(iC).psd(2,lf:rf));
            if od.S1.cell_type(iC) == 1
                msn_near_psd_difs = [msn_near_psd_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            else
                fsi_near_psd_difs = [fsi_near_psd_difs (od.S1.freqs(lf+p2)-od.S1.freqs(lf+p1))];
            end
        end
        
        lf = find(od.S2.freqs >= 40, 1, 'first');
        rf = find(od.S2.freqs <= 80, 1, 'last');
        for iC  = 1:length(od.S2.cell_type)
            if throw_msn2 == length(msn2) || od.S2.spec_res(iC).scount_ptile(1) < 200 || od.S2.spec_res(iC).scount_ptile(2) < 200
                continue
            end
            [~,p1] = max(od.S2.spec_res(iC).psd(1,lf:rf));
            [~,p2] = max(od.S2.spec_res(iC).psd(2,lf:rf));
            if od.S2.cell_type(iC) == 1
                msn_away_psd_difs = [msn_away_psd_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            else
                fsi_away_psd_difs = [fsi_away_psd_difs (od.S2.freqs(lf+p2)-od.S2.freqs(lf+p1))];
            end
        end
    end
end
figure;
subplot(2,1,1)
histogram(msn_away_psd_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_away_psd_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Away trials average LFP PSD');
subplot(2,1,2)
histogram(msn_near_psd_difs,-37.5:5:37.5, 'FaceColor','red','FaceAlpha',0.4)
hold on
histogram(fsi_near_psd_difs,-37.5:5:37.5, 'FaceColor','green','FaceAlpha',0.4)
title('Near trials average LFP PSD');
