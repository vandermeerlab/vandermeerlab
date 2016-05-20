function [keep_cycle, keep_idx] = AMPX_extract_max_cycle(cfg, data_filt)
%% AMPX_extract_max_cycle: find the maximum amplitude cycle in a filtered
% data event.  Finds the highest peak in the filtered data, then it finds
% the two adjacent troughs to complete one cycle around the center of the
% highest amplitude cycle.

% inputs:
%     - cfg [struct]: only needs the "cfg.chan" field
%     - filt_data [struct]: filtered data in the TSD format
%
%
% outputs:
%     - cycle_data [nChan x n]: matrix with rows for channels and all data
%       points within the cycle.
%
%
% EC 11/10/2015
%
% note: Since cycles are not always a fixed width across an event, this
% means that cycle_data can be different sizes across events.  Be prepared
% to deal with this.  I would recommend interolating to match the longest
% event.
%
% To Do: Could use a function that finds the best channel based on the BP power
% or the PSD for the band of interest if unfiltered.
%% extract the phase from the best channel

phase = -angle(hilbert(data_filt.data(cfg.chan,:)));

rmpath('/Users/ecarmichael/Documents/Github/fieldtrip/external/stats/vline.m')

%% find the peak in the power of the filtered singal.  Then get 1/2 a phase cycle on either side.

[pks, loc] = findpeaks(data_filt.data(cfg.chan, :), 'MinPeakDistance', cfg.MPD);%, 'MaxPeakWidth', cfg.MPW);


    peak = loc(pks==max(pks)); % get the max peak idx


if peak == loc(1)  % just in case it picks up the first peak which means it can't get the previous 180? point. 
    disp('First peak was the max, using second largest')
    ii = sort(pks);
    out = pks(pks == ii(end-1));

    peak = loc(pks == out);
end
phase_at_peak = phase(peak);
[~, xloc] = findpeaks(-1*data_filt.data(cfg.chan, :),'MinPeakDistance', cfg.MPD);
% next_peak = xloc(find(loc == peak)+1);
% prev_peak = xloc(find(loc == peak)+0);
[~, n_loc] = max(xloc>peak);
next_peak = xloc(n_loc);

[~, p_loc] = max(xloc(xloc<peak));
prev_peak = xloc(p_loc);

phase_at_next = phase( next_peak);
phase_at_prev = phase( prev_peak);

if isfield(cfg, 'plot') && strcmp(cfg.plot, 'on')
    h_cycle = figure;
    time_of_peak = data_filt.tvec(peak);
    subplot(411)
    hold on
    plot(data_filt.tvec, data_filt.data(cfg.chan, :))
    vline([data_filt.tvec(peak),data_filt.tvec(next_peak), data_filt.tvec(prev_peak)], {'r', 'b', 'g'}, {'peak', 'next', 'prev'})
    
    subplot(412)
    hold on
    plot(data_filt.tvec, data_filt.data(cfg.chan, :))
    vline([data_filt.tvec(peak),data_filt.tvec(next_peak), data_filt.tvec(prev_peak)], {'r', 'b', 'g'}, {'peak', 'next', 'prev'})
    xlim([data_filt.tvec(prev_peak)-(data_filt.tvec(next_peak) - data_filt.tvec(prev_peak)) data_filt.tvec(next_peak)+(data_filt.tvec(next_peak) - data_filt.tvec(prev_peak))]);
    
    subplot(413)
    hold on
    plot(data_filt.tvec, phase(:))
    vline([data_filt.tvec(peak),data_filt.tvec(next_peak), data_filt.tvec(prev_peak)], {'r', 'b', 'g'}, {'peak', 'next', 'prev'})
    xlim([data_filt.tvec(prev_peak)-(data_filt.tvec(next_peak) - data_filt.tvec(prev_peak)) data_filt.tvec(next_peak)+(data_filt.tvec(next_peak) - data_filt.tvec(prev_peak))]);
end

%% display the times
fprintf(['Time of peak:' num2str(data_filt.tvec(peak)) ' \nlength of cycle:' num2str(data_filt.tvec(next_peak)- data_filt.tvec(prev_peak)) '\n'])
fprintf(['Phase at start: ' num2str(phase_at_prev) '  Phase at peak:' num2str(phase_at_peak) '  Phase at end: ' num2str(phase_at_next) '\n'])

%% collect all the cycles across channels
for iChan = size(data_filt.data,1):-1:1
    keep_cycle(iChan, :) = data_filt.data(iChan, prev_peak:next_peak);
    keep_idx(iChan, :) = [prev_peak , next_peak];
    if isfield(cfg, 'plot') && strcmp(cfg.plot, 'on')
        subplot(414)
        hold on
        plot(data_filt.tvec(prev_peak: next_peak), keep_cycle(iChan,:))
        xlim([data_filt.tvec(prev_peak), data_filt.tvec(next_peak)])
    end
end

