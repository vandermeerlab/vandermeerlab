function evt = MatchGammaEvents(cfg_in,csc,lg_in,hg_in)
% function evt = MatchGammaEvents(cfg_in,csc,lg_in,hg_in)
%
% finds event times with low gamma power matching gamma events
%
% algorithm:
% for each incoming gamma event:
% - within window cfg.twin from incoming event center
% - find time with lowest gamma power that isn't already taken (see below)
% - when event found, apply mask so that subsequent events can't overlap
%   with already chosen events
%
% INPUT:
%
% csc: tsd with raw LFP to use for obtaining wideband gamma power
% lg_in, hg_in: times (in s), NOT iv's, of detected gamma events
%  NOTE cfg.evt_twin specifies width of events
%
% OUTPUT:
% evt.lg: center times for matched low gamma events
% evt.lg_dt: time differences from original to matched events (e.g. 2 means
%   matched event occurred 2 seconds after original event)
% evt.hg: center times for matched high gamma events
% evt.hg_dt: time differences from original to matched events (e.g. 2 means
%   matched event occurred 2 seconds after original event)

cfg_def = [];
cfg_def.twin = [-2 2]; % in s
cfg_def.f = [40 85]; % freq band for obtaining gamma power
cfg_def.evt_twin = [-0.2 0.2]; % size of event
cfg_def.verbose = 1;

cfg = ProcessConfig(cfg_def,cfg_in);

% first, filter LFP
cfg_filter = [];
cfg_filter.f = cfg.f;
cfg_filter.type = 'cheby1'; cfg_filter.order = 5;
gPow = FilterLFP(cfg_filter,csc);

% create running average of gamma amplitude, over window size corresponding
% to event (this way we can just find the min of the smoothed signal later)
gPow.data = abs(hilbert(gPow.data));
dt = median(diff(gPow.tvec));
span = round(sum(abs(cfg.evt_twin))./dt);
gPow.data = Match_filt_fnc(gPow.data, span);
gPow.data = smooth(gPow.data,span);

% initialize mask: can't use times that overlap with input events
mask = ones(size(gPow.data));
mask_tsd = tsd(gPow.tvec,mask);

% loop through events IN RANDOM ORDER! remember to return times
nlg = length(lg_in);
nhg = length(hg_in);
nall = nlg + nhg;

evt_t = cat(1,lg_in,hg_in); % times of original (input) events
evt_label = cat(1,zeros(nlg,1),ones(nhg,1)); % labels of original (input) events (0: lg, 1: hg)
evt_idx = cat(2,1:nlg,1:nhg); % idxs of original events
evt_matched = nan(size(evt_t)); % output times go here

evt_order = randperm(nall); % order to match events in.. loop through this from start to end

% mask out input events
for iEvt = 1:nall
   
    csc_idx = TSD_getidx2(csc,evt_t(iEvt)+2*cfg.evt_twin(1),evt_t(iEvt)+2*cfg.evt_twin(2));
    mask_tsd.data(csc_idx) = NaN;
    
end

% loop over input events to find matching time
for iEvt = 1:nall
    
   % pick up first event
   this_evt_t = evt_t(evt_order(iEvt));
    
   % for this event, grab piece of eligible filtered data
   this_tstart = this_evt_t + cfg.twin(1);
   this_tend = this_evt_t + cfg.twin(2);
   
   csc_idx = TSD_getidx2(csc,this_tstart,this_tend);
   
   this_gPow = gPow.data(csc_idx);
   this_mask = mask_tsd.data(csc_idx);
   this_tvec = gPow.tvec(csc_idx);
   
   this_gPow = this_gPow.*this_mask; % this NaNs out any unavailable time points
   
   % find minimum
   [~,min_idx] = min(this_gPow);
   matched_evt_t = this_tvec(min_idx);
   
   % update mask
   mask_idx = TSD_getidx2(csc,matched_evt_t+2*cfg.evt_twin(1),matched_evt_t+2*cfg.evt_twin(2));
   mask_tsd.data(mask_idx) = NaN;
   
   % figure out which event this was, anyway
   if evt_label(evt_order(iEvt)) == 0 % lg
       evt.lg(evt_idx(evt_order(iEvt))) = matched_evt_t;
       evt.lg_dt(evt_idx(evt_order(iEvt))) = matched_evt_t - this_evt_t;
   else % hg
       evt.hg(evt_idx(evt_order(iEvt))) = matched_evt_t;
       evt.hg_dt(evt_idx(evt_order(iEvt))) = matched_evt_t - this_evt_t;
   end
    
end