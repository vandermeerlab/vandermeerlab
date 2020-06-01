cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/psd_v1');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString1 = strcat(curRat,'*psd.mat');
    searchString2 = strcat(curRat,'*scount.mat');
    ofiles1 = dir(searchString1);
    ofiles2 = dir(searchString2);
    for jdx = 1:length(ofiles1)
        load(ofiles1(jdx).name);
        temp = load(ofiles2(jdx).name);
        for iC = 1:length(od.S1.cell_type)
            od.S1.spec_res(iC).scount_ptile = temp.od.S1.spec_res(iC).scount_ptile;
        end
        for iC = 1:length(od.S2.cell_type)
            od.S2.spec_res(iC).scount_ptile = temp.od.S2.spec_res(iC).scount_ptile;
        end
        save(ofiles1(jdx).name,'od'); % should add option to save in specified output dir
    end
end