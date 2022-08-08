function events_ts = LoadEvents(cfg_in)
% events = LoadEvents(cfg)
%
% generic loader for Neuralynx .nev files
%
% INPUTS:
%
% cfg.eventList: cell array with event strings to process
%  if not specified, load all
% cfg.eventLabel: cell array with strings to use as labels
%  (should match eventList in size)
%
% OUTPUTS:
%
% events_ts: ts object with event timestamps
%
% MvdM 2014-06-18, 25 (use cfg_in)

cfg_def = [];
cfg = ProcessConfig(cfg_def,cfg_in); % this takes fields from cfg_in and puts them into cfg

mfun = mfilename;

% load events file
fn = FindFile('*Events.nev');
if isempty(fn)
   error('LoadEvents: no events file found.'); 
end

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

% obtain event labels to process
if ~isfield(cfg,'eventList')
        
        cfg.eventList = unique(EventStrings);
        
else
        
        if ~isa(cfg.eventList,'cell')
            error('LoadEvents: cfg.eventList should be a cell array.');
        end
        
end

if ~isfield(cfg,'eventLabel') % if no label defined, just use the list 
        
        cfg.eventLabel = cfg.eventList;
        
end

% get times for each event label
events_ts = ts;

for iEvent = 1:length(cfg.eventList)
   
    ev_string = cfg.eventList{iEvent};

    ev_id = strncmp(ev_string,EventStrings,length(ev_string));
    ev_t = EVTimeStamps(ev_id)*10^-6;
    
    % check if this eventLabel already exists, if so append (not create new)
    label_idx = strmatch(cfg.eventLabel{iEvent},events_ts.label,'exact');
    
    if isempty(label_idx) % label doesn't exist, create new
        
        events_ts.t{iEvent} = ev_t;
        events_ts.label{iEvent} = cfg.eventLabel{iEvent};
        
    else % label exists, append and sort (sort is helpful for later use with nearest_idx3)
        
        events_ts.t{label_idx} = sort(cat(2,events_ts.t{label_idx},ev_t));
                
    end
    
end

% check if ExpKeys available
keys_f = FindFiles('*keys.m');
if ~isempty(keys_f)
    run(keys_f{1});
    events_ts.cfg.ExpKeys = ExpKeys;
end

% add sessionID
[~,events_ts.cfg.SessionID,~] = fileparts(pwd);

events_ts.cfg.filename = fn;

% housekeeping
events_ts.cfg.history.mfun = cat(1,events_ts.cfg.history.mfun,mfun);
events_ts.cfg.history.cfg = cat(1,events_ts.cfg.history.cfg,{cfg});
