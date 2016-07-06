function [CSD] = AMPX_cycle_kCSD(cfg_in, all_cycles)
%% AMPX_cycle_kCSD: computes the kernal current source density (kCSD) for
%  each three cycle event in all_cycles.
%
%   *kCSD toolbox required: http://software.incf.org/software/kcsd
%
%          Inputs:
%           - cfg_in: [struct] contains paramters
%           - all_cycles: [struct] output from AMPX_extract_max_cycle
%          Outputs:
%           - cycle_csd: [struct] contains the
%
% EC - 2016-05-19


%% parameters
% cfg_in.index = []; % which events to use.  If empty then all will be used.
cfg_def.buffer = 0.2; % time (s) to include before and after gamma
cfg_def.interp_factor = 0.1; % new spacing for interpolation
cfg_def.chanlabel = diag(reshape(1:64,8,8)); % vertical section with no bad channels to plot for raw data.  Default for A8x8 AMPX
cfg_def.scalefactor = 600; % scale factor for raw data plot
cfg_def.Fsvid = 10; % sampling frequency for making the video
cfg_def.clims = [-800 1000]; % range of colour scale in LFP image
cfg_def.clims_csd = [-700 700]; % range of colour scale in CSD image
cfg_def.font_size = 13; % fontsize
cfg_def.t_font_size = 14; % fontsize for title
cfg_def.ref_chan = 57;
cfg_def.data_smooth = 10;
cfg_def.session_type = 'pre';
cfg_def.smooth = 10;
cfg_def.plot = 'on';
cfg_def.peak_to_use = 'next';
cfg = ProcessConfig2(cfg_def,cfg_in);

%% cycle all the event
peaks = {'prev', 'center', 'next'};

for ievt =length(all_cycles.peaks.center.data):-1:1;
    for iPeak = 1:length(peaks)
        cfg_kCSD = [];
        cfg_kCSD.R = 0.2;
        csd_struct = AMPX_kCSD(cfg_kCSD, all_cycles.peaks.(peaks{iPeak}).tvec{ievt}, all_cycles.peaks.(peaks{iPeak}).data{ievt}, all_cycles.ExpKeys);
        CSD.(peaks{iPeak}).all_csd_1d(:,:,ievt) = csd_struct.csd_1d;
        CSD.(peaks{iPeak}).cfg{ievt} = csd_struct.cfg;
        CSD.(peaks{iPeak}).R(ievt) = csd_struct.cfg.R;
        axis off
    end
end
%     if CSD.(cfg.band).failures(1, ievt) ==0
%         mkdir([datestr(date, 'yyyy-mm-dd') '/CSD/'])
%         saveas(h_cycle, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_' num2str(ievt) '_CSD_cycle.fig'])
%         saveas(h_cycle, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_' num2str(ievt) '_CSD_cycle.png'])
%     end
close all


%% plot the average on a new figure and save the data.
h_avg = figure(10);
h1 = subplot(1,3,1);
set(gcf, 'position', [677 526 736 215])
imagesc(mean(CSD.prev.all_csd_1d, 3))
p = get(h1, 'pos');
p(3) = p(3) + 0.05;
set(h1, 'pos', p);
% text(30, 5, 'prev')
axis off
h2 = subplot(1,3,2);
imagesc(mean(CSD.center.all_csd_1d, 3))
% text(30, 5, 'center')
axis off
p = get(h2, 'pos');
p(3) = p(3) + 0.05;
set(h2, 'pos', p);


h3 = subplot(1,3,3);
imagesc(mean(CSD.next.all_csd_1d, 3))
% text(30, 5, 'next')
axis off
p = get(h3, 'pos');
p(3) = p(3) + 0.05;
set(h3, 'pos', p);


%%
print(['C:\Users\mvdmlab\Dropbox\Naris_Paper\CSD_figure\CSD_' cfg.band 'cycles'], -depsc)
saveas(h_avg, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_avg_CSD_cycle.fig'])
saveas(h_avg, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_avg_cycle.png'])
saveas(h_avg, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_avg_cycle.png'])
