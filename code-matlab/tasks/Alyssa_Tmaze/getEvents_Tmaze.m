function events = getEvents_Tmaze(varargin)
% events = getEvents_Tmaze(varargin)
%
% get events for Alyssa's T-maze task
%
% needs to be run from data folder! 
%
% varargins:
%
% rawTimeStamps = 0;
%
% MvdM 2013-08-16 initial version


rawTimeStamps = 0;
extract_varargin;

eventList = {'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0020).','center_pb',...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0080).','left_pb',...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).','right_pb',...
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0004).','food_dispensed',...
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0040).','water_dispensed'};


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