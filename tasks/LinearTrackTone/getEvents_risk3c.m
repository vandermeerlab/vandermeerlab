function events = getEvents_risk3c(varargin)
% events = getEvents_risk3c(varargin)
%
% get events for LinearTrackTone task (risk version, 3 cues)
%
% MvdM 06-13-2012

eventList = {'1 or 5 pellet cue','hirisk',...
    '2 or 4 pellet cue','lorisk',...
    '3 pellet cue','norisk',...
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
