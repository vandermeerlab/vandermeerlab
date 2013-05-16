function [trl, event] = ft_trialfun_lineartracktone(cfg)
% function [trl, event] = ft_trialfun_lineartracktone(cfg)
%
% basic trialfun for lineartracktone
%
% needs cfg.trialdef.hdr in order to convert to sample idx
% needs cfg.trialdef.eventtype to determine what to extract
% needs cfg.trialdef.pre and cfg.trialdef.post
% can use cfg.trialdef.eventvalue % dependent on eventtype
% can use cfg.trialdef.block % {'value','risk','both'}
%
% MvdM 2013-04-20

evt = getEvents_LinearTrackTone('rawTimeStamps',1); % get all events

switch cfg.trialdef.eventtype
    
    case 'nosepoke'
        
        event.ts = cat(2,evt.n0,evt.n1);
        
        %if cfg.trialdef.eventvalue == 0 % feeder 0
            
        %else % feeder 1
            
        %end
    
end

% restrict to block of interest, if applicable


event.ts = sort(event.ts,'ascend');

% convert to sample idx based on information in header, i.e. in nlx
% timestamp space
tvec = cfg.trialdef.hdr.FirstTimeStamp:cfg.trialdef.hdr.TimeStampPerSample:cfg.trialdef.hdr.LastTimeStamp;
tvec = double(tvec);

for iE = length(event.ts):-1:1
    
   event.idx(iE) = nearest(tvec,event.ts(iE));
    
end

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * cfg.trialdef.hdr.Fs);
posttrig =  round(cfg.trialdef.post * cfg.trialdef.hdr.Fs);

% construct trl output
trl = [];
for j = 1:numel(event.idx)

    % find the sample corrsponding to this timestamp
    
    trlbegin = event.idx(j) + pretrig;       
    trlend   = event.idx(j) + posttrig;       
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];

end