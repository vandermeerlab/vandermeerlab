%% Script to generate plots
% 
% clear;
% close all;
% load('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/sts_v4/R117-2007-06-01_sts.mat');
% c_list = {'k','b','g','r'};
% % for near-trials;
% for iC  = 1:length(od.S1.cell_type)
%     close;
%     legend_text = cell(1,4);
%     for iP = 1:4
%         subplot(2,1,1);
%         plot(od.S1.freqs,10*log(od.S1.spec_res(iC).sta_spec_ptile(iP,:)), c_list{iP});
%         hold on;
%         plot(od.S1.freqs,10*log(od.S1.spec_res(iC).w_sta_spec_ptile(iP,:)),cat(2,':',c_list{iP}));
%         hold on;
%         plot(od.S1.freqs,10*log(od.S1.spec_res(iC).uw_sta_spec_ptile(iP,:)),cat(2,'--',c_list{iP}));
%         hold on;
%         subplot(2,1,2);
%         plot(od.S1.freqs,10*log(od.S1.spec_res(iC).sts_ptile(iP,:)),c_list{iP});
%         hold on;
%         legend_text{iP} = cat(2, num2str(od.S1.spec_res(iC).ptile_mfrs(iP,1)), ...
%             ' Hz - ', num2str(od.S1.spec_res(iC).ptile_mfrs(iP,2)), ' Hz, Spikes: ', num2str(od.S1.spec_res(iC).scount_ptile(iP)));
%         dummy = 1;
%     end
%     subplot(2,1,1);
%     title('STA spectrum for away trials','FontSize', 18);
%     subplot(2,1,2);
%     title('Average STS for away trials','FontSize', 18);
%     legend(legend_text, 'FontSize', 14);
%     dummy = 2;   
% end



cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/sts_v4');
rats = {'R117'};%,'R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*sts.mat');
    ofiles = dir(searchString);
    c_list = {'b','r'};
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        
        % for near-trials;
        for iC  = 1:length(od.S1.cell_type)
            legend_text = cell(1,2);
            o_prefix = extractBefore(od.S1.label{iC},'.t');
            fig = figure;
            for iP = 1:2
                subplot(3,2,iP);
                wl = floor(length(od.S1.spec_res(1).sta_ptile(iP,:))/2);
                plot((-wl:1:wl),od.S1.spec_res(1).sta_ptile(iP,:),c_list{iP});
                xlim([-wl wl]);
                subplot(3,2,3:4);
                plot(od.S1.freqs,10*log(od.S1.spec_res(iC).sta_mtspec_ptile(iP,:)),c_list{iP});
                hold on;
                subplot(3,2,5:6);
                plot(od.S1.freqs,10*log(od.S1.spec_res(iC).mtsts_ptile(iP,:)),c_list{iP});
                hold on;
                legend_text{iP} = cat(2, num2str(od.S1.spec_res(iC).ptile_mfrs(iP,1)), ...
                    ' Hz - ', num2str(od.S1.spec_res(iC).ptile_mfrs(iP,2)), ' Hz, Spikes: ', ...
                    num2str(od.S1.spec_res(iC).scount_ptile(iP)));
            end
            subplot(3,2,1);
            title('STA for Low FR near trials','FontSize', 12);
            subplot(3,2,2);
            title('STA for High FR near trials','FontSize', 12);
            subplot(3,2,3:4);
            title('STA spectrum for near trials','FontSize', 12);
            legend(legend_text);
            subplot(3,2,5:6);
            title('Average STS for near trials','FontSize', 12);
            legend(legend_text);
            if od.S1.cell_type(iC) == 1
                o_name = cat(2, o_prefix,'_near_MSN');
            else
                o_name = cat(2, o_prefix,'_near_FSI');
            end
            WriteFig(fig,o_name,1);
            close all; 
        end
        
        % for away-trials;
        for iC  = 1:length(od.S2.cell_type)
            legend_text = cell(1,2);
            o_prefix = extractBefore(od.S2.label{iC},'.t');
            fig = figure;
            for iP = 1:2
                subplot(3,2,iP);
                wl = floor(length(od.S2.spec_res(1).sta_ptile(iP,:))/2);
                plot((-wl:1:wl),od.S2.spec_res(1).sta_ptile(iP,:),c_list{iP});
                xlim([-wl wl]);
                subplot(3,2,3:4);
                plot(od.S2.freqs,10*log(od.S2.spec_res(iC).sta_mtspec_ptile(iP,:)),c_list{iP});
                hold on;
                subplot(3,2,5:6);
                plot(od.S2.freqs,10*log(od.S2.spec_res(iC).mtsts_ptile(iP,:)),c_list{iP});
                hold on;
                legend_text{iP} = cat(2, num2str(od.S2.spec_res(iC).ptile_mfrs(iP,1)), ...
                    ' Hz - ', num2str(od.S2.spec_res(iC).ptile_mfrs(iP,2)), ' Hz, Spikes: ', ...
                    num2str(od.S2.spec_res(iC).scount_ptile(iP)));
            end
            subplot(3,2,1);
            title('STA for Low FR away trials','FontSize', 12);
            subplot(3,2,2);
            title('STA for High FR away trials','FontSize', 12);
            subplot(3,2,3:4);
            title('STA spectrum for away trials','FontSize', 12);
            legend(legend_text);
            subplot(3,2,5:6);
            title('Average STS for away trials','FontSize', 12);
            legend(legend_text);
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



        
