function events = getEvents_LinearTrackTone(varargin)
% events = getEvents_LinearTrackTone(varargin)
%
% get events for LinearTrackTone task
%
% needs to be run from data folder! 
%
% varargins:
%
% eventList = {'1 pellet cue','c1',... where first element is the
% EventString and second is the field resulting time(stamp)s are stored
%
% rawTimeStamps = 0;
%
% MvdM 05-05-2012 initial version
% MvdM 04-20-2013 edit to interface with ft trialfun

eventList = {'1 pellet cue','c1',...
    '3 pellet cue','c3',...
    '5 pellet cue','c5',...
    '1 pellet dispensed','d1',...
    '3 pellet dispensed','d3',...
    '5 pellet dispensed','d5',...
    'Feeder 0 nosepoke','n0',...
    'Feeder 1 nosepoke','n1'};

rawTimeStamps = 0;
extract_varargin;

fn = FindFile('*Events.nev');

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

events = [];

for iEvent = 1:2:length(eventList)
   
    ev_string = eventList{iEvent};
    ev_target = eventList{iEvent+1};
    
    ev_id = strncmp(ev_string,EventStrings,length(ev_string));
    
    if rawTimeStamps
        ev_t = EVTimeStamps(ev_id);
    else
        ev_t = EVTimeStamps(ev_id)*10^-6;
    end
    
    events.(ev_target) = ev_t;
    
end