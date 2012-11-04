function events = getEvents(varargin)
% events = getEvents(varargin)
%
% get events for LinearTrackTone task
%
% MvdM 05-05-2012

eventList = {'1 pellet cue','c1',...
    '3 pellet cue','c3',...
    '5 pellet cue','c5',...
    '1 pellet dispensed','d1',...
    '3 pellet dispensed','d3',...
    '5 pellet dispensed','d5',...
    'Feeder 0 nosepoke','n0',...
    'Feeder 1 nosepoke','n1'};

eventListAlt = {0,'n0',1,'n1',20,'c1',21,'c3',22,'c5',30,'d1',31,'d3',32,'d5'}; % not used, for checking...
    
getPredictionErrorEvents = 0;
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

eventsAlt = [];

for iEvent = 1:2:length(eventList)
   
    ev_string = eventListAlt{iEvent};
    ev_target = eventListAlt{iEvent+1};
    
    ev_id = find(EventIDs == ev_string);
    ev_t = EVTimeStamps(ev_id)*10^-6;
    
    eventsAlt.(ev_target) = ev_t;
    
end