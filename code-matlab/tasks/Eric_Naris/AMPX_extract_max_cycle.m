function [cycle, h_cycle] = AMPX_extract_max_cycle(cfg, data_filt)
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
use_flag = 0; % used to detect events with abnormal shape.
phase = -angle(hilbert(data_filt.data(cfg.chan,:)));

rmpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip\external\stats\')

%% find the peak in the power of the filtered singal.  Then get 1/2 a phase cycle on either side.

[pks, loc] = findpeaks(data_filt.data(cfg.chan, :), 'MinPeakDistance', cfg.MPD);%, 'MaxPeakWidth', cfg.MPW);

[pk_in, loc_in] = findpeaks(-data_filt.data(cfg.chan, :), 'MinPeakDistance', cfg.MPD, 'MinPeakHeight', max(pks)*0.8);

    peak = loc(pks==max(pks)); % get the max peak idx

%%
if peak < (length(data_filt.data(cfg.chan, :))/20) || peak > (length(data_filt.data(cfg.chan, :)) - (length(data_filt.data(cfg.chan, :))/20))%|| loc(1)<loc_in(1)) % just in case it picks up the first peak which means it can't get the previous 180? point.
    disp('Largest peak was in the buffer')
    ii = sort(pks);
    out = pks(pks == ii(end-1));
    use_flag = 1;
    peak = loc(pks == out);
    cycle = [];
    h_cycle = [];
    return
end
phase_at_peak = phase(peak);
[~, xloc] = findpeaks(-1*data_filt.data(cfg.chan, :),'MinPeakDistance', cfg.MPD);
% next_peak = xloc(find(loc == peak)+1);
% prev_peak = xloc(find(loc == peak)+0);
[~, n_loc] = max(xloc>peak);
next_val = xloc(n_loc);

[~, p_loc] = max(xloc(xloc<peak));
prev_val = xloc(p_loc);

phase_at_next = phase( next_val);
phase_at_prev = phase( prev_val);

% find the valley before the prev peak and after the next peak. 
[~, ab_loc] = findpeaks(abs(data_filt.data(cfg.chan, :)),'MinPeakDistance', round(cfg.MPD/2));

% get the previous cycle
prev_loc = find(loc ==peak);
prev_peak = loc(prev_loc -1);
% prev_peak_val = xloc(prev_loc -2);
prev_loc_ab = find(ab_loc == prev_peak);
if prev_loc_ab == 1
        cycle = [];
    h_cycle = [];
    disp('####################Event does not meet criteria#######################')
    return
end
prev_peak_val = ab_loc(prev_loc_ab-1);

% get the next cycle
next_loc = find(loc ==peak);
if next_loc == length(ab_loc)
        cycle = [];
    h_cycle = [];
    disp('####################Event does not meet criteria#######################')
    return
end
next_peak = loc(next_loc +1);
next_loc_ab = find(ab_loc == next_peak);
if next_loc_ab == length(ab_loc)
        cycle = [];
    h_cycle = [];
    disp('####################Event does not meet criteria#######################')
    return
end
next_peak_val = ab_loc(next_loc_ab+1);


%% check if there is an issue with the event (not a previous or next peak
% within the event.
if prev_peak_val > prev_peak || next_peak_val < next_peak
    cycle = [];
    h_cycle = [];
    disp('####################Event does not meet criteria#######################')
    return
end

%% plot

% if isfield(cfg, 'plot') && strcmp(cfg.plot, 'on')
    %%
    h_cycle = figure;
    time_of_peak = data_filt.tvec(peak);
    subplot(5,3,1:3)
    hold on
    plot(data_filt.tvec, data_filt.data(cfg.chan, :))
    vline([data_filt.tvec(peak),data_filt.tvec(next_val), data_filt.tvec(prev_val)], {'-r', '-b', '-g'}, {'peak', 'next', 'prev'})
    vline([data_filt.tvec(prev_peak), data_filt.tvec(prev_peak_val)], {'--c', '--k'}, {'   peak_peak', '   prev_val'})
    vline([data_filt.tvec(next_peak), data_filt.tvec(next_peak_val)], {'--m', '--g'}, {'   next_peak', '   next_val'})


    subplot(5,3,4:6)
    hold on
    plot(data_filt.tvec, data_filt.data(cfg.chan, :))
    vline([data_filt.tvec(peak),data_filt.tvec(next_val), data_filt.tvec(prev_val)], {'-r', '-b', '-g'}, {'peak', 'next', 'prev'})
    vline([data_filt.tvec(prev_peak), data_filt.tvec(prev_peak_val)], {'--c', '--k'}, {'   peak_peak', '   prev_val'})
    vline([data_filt.tvec(next_peak), data_filt.tvec(next_peak_val)], {'--m', '--g'}, {'   next_peak', '   next_val'})
    xlim([data_filt.tvec(prev_val)-(data_filt.tvec(next_val) - data_filt.tvec(prev_val)) data_filt.tvec(next_val)+(data_filt.tvec(next_val) - data_filt.tvec(prev_val))]);
    
    subplot(5,3,7:9)
    hold on
    plot(data_filt.tvec, phase(:))
    vline([data_filt.tvec(peak),data_filt.tvec(next_val), data_filt.tvec(prev_val)], {'r', 'b', 'g'}, {'peak', 'next', 'prev'})
    xlim([data_filt.tvec(prev_val)-(data_filt.tvec(next_val) - data_filt.tvec(prev_val)) data_filt.tvec(next_val)+(data_filt.tvec(next_val) - data_filt.tvec(prev_val))]);
% end

%% display the times
fprintf(['Time of peak:' num2str(data_filt.tvec(peak)) ' \nlength of cycle:' num2str(data_filt.tvec(next_val)- data_filt.tvec(prev_val)) '\n'])
fprintf(['Phase at start: ' num2str(phase_at_prev) '  Phase at peak:' num2str(phase_at_peak) '  Phase at end: ' num2str(phase_at_next) '\n'])

%% collect all the cycles across channels
keep_idx = [prev_val , next_val];
keep_idx_prev = [prev_peak_val, prev_val];
keep_idx_next = [next_val, next_peak_val];

for iChan = size(data_filt.data,1):-1:1
    keep_cycle(iChan, :) = data_filt.data(iChan, prev_val:next_val);
    if isfield(cfg, 'plot') && strcmp(cfg.plot, 'on')
        subplot(5,3,11)
        hold on
        plot(data_filt.tvec(prev_val:next_val), keep_cycle(iChan,:), 'r')
        xlim([data_filt.tvec(prev_val) data_filt.tvec(next_val)])
    end
end
xlabel('Center')

% collect the prev_peak
for iChan = size(data_filt.data,1):-1:1
    keep_cycle_prev(iChan, :) = data_filt.data(iChan, prev_peak_val:prev_val);
    if isfield(cfg, 'plot') && strcmp(cfg.plot, 'on')
        subplot(5,3,10)
        hold on
        plot(data_filt.tvec(prev_peak_val:prev_val), keep_cycle_prev(iChan,:), 'b')
        xlim([data_filt.tvec(prev_peak_val) data_filt.tvec(prev_val)])
    end
end
xlabel('Prev')


% collect the next_peak
for iChan = size(data_filt.data,1):-1:1
    keep_cycle_next(iChan, :) = data_filt.data(iChan, next_val:next_peak_val);
    if isfield(cfg, 'plot') && strcmp(cfg.plot, 'on')
        subplot(5,3,12)
        hold on
        plot(data_filt.tvec(next_val: next_peak_val), keep_cycle_next(iChan,:), 'g')
        xlim([data_filt.tvec(next_val) data_filt.tvec(next_peak_val)])
    end
end
xlabel('Next')
%% gather the outputs
cycle.center.cycle = keep_cycle;
cycle.center.idx = keep_idx;

cycle.prev.cycle = keep_cycle_prev;
cycle.prev.idx = keep_idx_prev;

cycle.next.cycle = keep_cycle_next;
cycle.next.idx = keep_idx_next;

%% check for the normal shape of the event.  Needs to have a center and a
% peak on either side.
peak_idx= find(loc == peak);
if peak_idx <2
    use_flag = 2;
elseif peak_idx <2 && use_flag ~=0
    use_flag = 3;
end
addpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip\external\stats\')

