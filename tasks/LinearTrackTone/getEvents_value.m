function events = getEvents_value(varargin)
% events = getEvents_value(varargin)
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

convertToSamples = 0;
getPredictionErrorEvents = 0;
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
    elseif convertToSamples
        ev_t = EVTimeStamps(ev_id)*10^-6;
        for j = numel(ev_t):-1:1
            ev_t(j) = nearest(tvec,ev_t(j));
        end
    else
        ev_t = EVTimeStamps(ev_id)*10^-6;
    end
    
    events.(ev_target) = ev_t;
    
end