function events = getEvents_risk(varargin)
% events = getEvents_risk(varargin)
%
% get events for LinearTrackTone task (risk version, 2 cues)
%
% MvdM 05-15-2012

eventList = {'1 or 5 pellet cue','risky',...
    '3 pellet cue','nonrisky',...
    '1 pellet dispensed','d1',...
    '3 pellet dispensed','d3',...
    '5 pellet dispensed','d5',...
    'Feeder 0 nosepoke','n0',...
    'Feeder 1 nosepoke','n1'};

extract_varargin;

fn = FindFile('*Events.nev');

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

events = [];

for iEvent = 1:2:length(eventList)
   
    ev_string = eventList{iEvent};
    ev_target = eventList{iEvent+1};
    
    ev_id = strncmp(ev_string,EventStrings,length(ev_string));
    ev_t = EVTimeStamps(ev_id)*10^-6;
    
    events.(ev_target) = ev_t;
    
end
