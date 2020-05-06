%% Script to generate plots

clear;
close all;
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp/STS_results');
rats = {'R117','R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        % generate plots for near trials
        for iC = 1:length(od.S1.t)
            o_prefix = extractBefore(od.S1.label{iC},'.t');
            fig = figure;
            subplot(2,1,1);
            for iT = 1:length(od.S1.trial_starts)
                plot(od.S1.freqs, 10*log(od.S1.sts{iC}(iT,:)), 'Color', 'blue')
                hold on
            end
            title('Average STS for near trials','FontSize', 18);
            subplot(2,1,2);
            histogram(od.S1.mfr{iC});
            title('Distribution of MFR over near trials','FontSize', 18);
            if od.S1.cell_type(iC) == 1
                o_name = cat(2, o_prefix,'_near_MSN');
            else
                o_name = cat(2, o_prefix,'_near_FSI');
            end
            WriteFig(fig,o_name,1);
            close all;
        end
        % generate plots for away trials
        for iC = 1:length(od.S2.t)
            o_prefix = extractBefore(od.S2.label{iC},'.t');
            fig = figure;
            subplot(2,1,1);
            for iT = 1:length(od.S2.trial_starts)
                plot(od.S2.freqs, 10*log(od.S2.sts{iC}(iT,:)), 'Color', 'blue')
                hold on
            end
            title('Average STS for away trials','FontSize', 18);
            subplot(2,1,2);
            histogram(od.S2.mfr{iC});
            title('Distribution of MFR over away  trials','FontSize', 18);
            if od.S2.cell_type(iC) == 1
                o_name = cat(2, o_prefix,'_away_MSN');
            else
                o_name = cat(2, o_prefix,'_away_FSI');
            end
            WriteFig(fig,o_name,1);
            close all;
        end
    end
end



        
