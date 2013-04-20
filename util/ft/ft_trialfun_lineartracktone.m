function [trl, event] = ft_trialfun_mvdmlab_nosepoke(cfg)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
%event = ft_read_event(cfg.dataset);

% search for "trigger" events
%value  = [event(find(strcmp('trigger', {event.type}))).value]';
%sample = [event(find(strcmp('trigger', {event.type}))).sample]';

% nosepokes (all)
events = getEvents_value('convertToSamples',1,'tvec',cfg.tvec,'eventList',{'Feeder 0 nosepoke','n0','Feeder 1 nosepoke','n1'});
evt_timestamps = cat(2,events.n0,events.n1); % nosepokes, both feeders

% cues (5)
events = getEvents_value('convertToSamples',1,'tvec',cfg.tvec,'eventList',{'5 pellet cue','c5'});
evt_timestamps = events.c5;

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
posttrig =  round(cfg.trialdef.post * hdr.Fs);

% look for the combination of a trigger "7" followed by a trigger "64" 
% for each trigger except the last one
trl = [];
for j = 1:numel(evt_timestamps)

    % find the sample corrsponding to this timestamp
    
    trlbegin = evt_timestamps(j) + pretrig;       
    trlend   = evt_timestamps(j) + posttrig;       
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];

end

event = events;