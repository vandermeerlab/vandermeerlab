cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString1 = strcat(curRat,'*full.mat');
    searchString2 = strcat(curRat,'*sts.mat');
    ofiles1 = dir(searchString1);
    ofiles2 = dir(searchString2);
    for jdx = 1:length(ofiles1)
        load(ofiles1(jdx).name);
        temp = load(ofiles2(jdx).name);
        msn = find(temp.od.S1.cell_type==1);
        fsi = find(temp.od.S1.cell_type==2);
        %For near trials
        for iC = 1:length(msn)
            od.S1.msn_res(iC).sta_full = temp.od.S1.spec_res(msn(iC)).sta_full;
            od.S1.msn_res(iC).mtsts_full = temp.od.S1.spec_res(msn(iC)).mtsts_full;
            od.S1.msn_res(iC).sta_mtspec_full = temp.od.S1.spec_res(msn(iC)).sta_mtspec_full;
%             od.S1.msn_res(iC).psd = temp.od.S1.spec_res(msn(iC)).psd;
        end
        for iC = 1:length(fsi)
            od.S1.fsi_res(iC).sta_full = temp.od.S1.spec_res(fsi(iC)).sta_full;
            od.S1.fsi_res(iC).mtsts_full = temp.od.S1.spec_res(fsi(iC)).mtsts_full;
            od.S1.fsi_res(iC).sta_mtspec_full = temp.od.S1.spec_res(fsi(iC)).sta_mtspec_full;
            od.S1.fsi_res(iC).unsampled_sta_ptile = temp.od.S1.spec_res(fsi(iC)).sta_ptile;
            od.S1.fsi_res(iC).unsampled_mtsts_ptile = temp.od.S1.spec_res(fsi(iC)).mtsts_ptile;
            od.S1.fsi_res(iC).unsampled_sta_mtspec_ptile = temp.od.S1.spec_res(fsi(iC)).sta_mtspec_ptile;
%             od.S1.fsi_res(iC).psd = temp.od.S1.spec_res(fsi(iC)).psd;
        end
        msn = find(temp.od.S2.cell_type==1);
        fsi = find(temp.od.S2.cell_type==2);
        %For away trials
        for iC = 1:length(msn)
            od.S2.msn_res(iC).sta_full = temp.od.S2.spec_res(msn(iC)).sta_full;
            od.S2.msn_res(iC).mtsts_full = temp.od.S2.spec_res(msn(iC)).mtsts_full;
            od.S2.msn_res(iC).sta_mtspec_full = temp.od.S2.spec_res(msn(iC)).sta_mtspec_full;
%              od.S2.msn_res(iC).psd = temp.od.S2.spec_res(msn(iC)).psd;
        end
        for iC = 1:length(fsi)
            od.S2.fsi_res(iC).sta_full = temp.od.S2.spec_res(fsi(iC)).sta_full;
            od.S2.fsi_res(iC).mtsts_full = temp.od.S2.spec_res(fsi(iC)).mtsts_full;
            od.S2.fsi_res(iC).sta_mtspec_full = temp.od.S2.spec_res(fsi(iC)).sta_mtspec_full;
            od.S2.fsi_res(iC).unsampled_sta_ptile = temp.od.S2.spec_res(fsi(iC)).sta_ptile;
            od.S2.fsi_res(iC).unsampled_mtsts_ptile = temp.od.S2.spec_res(fsi(iC)).mtsts_ptile;
            od.S2.fsi_res(iC).unsampled_sta_mtspec_ptile = temp.od.S2.spec_res(fsi(iC)).sta_mtspec_ptile;
%             od.S2.fsi_res(iC).psd = temp.od.S2.spec_res(fsi(iC)).psd;
        end
        
        save(ofiles1(jdx).name,'od'); % should add option to save in specified output dir
    end
end