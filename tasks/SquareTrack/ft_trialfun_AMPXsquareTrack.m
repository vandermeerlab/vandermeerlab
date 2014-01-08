function [trl, event] = ft_trialfun_AMPXsquareTrack(cfg)
% [trl, event] = ft_trialfun_AMPXsquareTrack(cfg)
%
% basic trialfun for Eric's square track
%
% uses previously saved *events/m file
%
% needs cfg.trialdef.hdr in order to convert to sample idx
% needs cfg.trialdef.eventtype to determine what to extract
% needs cfg.trialdef.pre and cfg.trialdef.post
% needs cfg.trialdef.eventtype {'fires','no_fires','nose_in_fires','nose_in_nofires'}
% can use cfg.trialdef.eventID [n x 1 double] with IDs of sites to include
% [1 2 3 4], default all
%
% MvdM 2013-12-13

load(FindFile('*events.mat'));

% generate tvec to convert times to idx
tvec = 0:1./cfg.trialdef.hdr.Fs:cfg.trialdef.hdr.filelength_sec;
tvec = tvec(1:end-1)';

event.ts = getfield(evt,cfg.trialdef.eventtype);

% if selected, only get certain events
if isfield(cfg.trialdef,'eventID')
    [~,idx,~] = intersect(event.ts.ID,cfg.trialdef.eventID);
    event.ts = event.ts.t(idx);
else % include all
    event.ts = event.ts.t;
end
    
event.ts = sort(event.ts,'ascend');

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