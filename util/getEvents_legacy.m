function events = getEvents_legacy(varargin)
% events = getEvents_legacy(varargin)
%
% get events for legacy Redish lab tasks
%
% MvdM 05-20-2012

extract_varargin;

eventList = {32768,'F1',32771,'F1',32769,'F2',32770,'F2'};
fn = FindFile('*Events.nev');

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

events = [];

for iEvent = 1:2:length(eventList)
   
    ev_string = eventList{iEvent};
    ev_target = eventList{iEvent+1};
    
    ev_id = find(TTLs == ev_string);
    ev_t = EVTimeStamps(ev_id)*10^-6;
    
    if ~isempty(ev_t)
        % remove doubles
        min_t = 0.5; % seconds
        keep = cat(2,1,(diff(ev_t) >= min_t));
        ev_t = ev_t(logical(keep));
        
        events.(ev_target) = ev_t;
    end
    
end
