function [all_events, all_events_str] = AMPX_format_task_events(data, evt, varargin)
%AMPX_format_task_events this function will take the data that has already
%been preprocessed and locate only the events of interest. The data will
%then be recompiled into a single data stream that will allow for gamma
%detection and analysis.  

% INPUTS:
% data [struct]: contains the hdr, tvec, channels, labels
% evt [strcut]: the event file that should have been made during task preprocessing. this will contain all the events for each feeder fire plus omissions (later it will include photobeams.)
% event_ID [array]: will contain the feeder IDs of interest (below)
%Site IDs: feeder 1 = 1; 
%          feeder 2 = 2; 
%          feeder 3 = 3; 
%          feeder 4 = 4; 
%          Omissions = 5; 
%
%OUTPUTS
% data [struct] contains all the same fields but the channels and tvec will
% be altered to only contain the events of interest as they have been
% recompiled. 

%% initialize
load('events.mat')
event_ID = 1:5;
extract_varargin

if isfield(evt, 'task') == 0
    task = AMPX_task_events(data);
    evt.task = task.task;
end
%% determine the events of interest using the event_ID
% if isempty(event_ID) ~=1
    all_times = []; all_IDs = []; all_events =[];
    for iIDs = event_ID
%         disp(iIDs)
        e = ['e' num2str(iIDs)];
        events_IDs = evt.task.ids == iIDs;
        events_times = evt.task.times(events_IDs);
        all_times = [all_times; events_times];
        all_events_str.(e) = events_times;
        all_IDs = [all_IDs; ones(length(events_times),1)*iIDs];
    end
    all_events(:,1) = all_times; 
    all_events(:,2) = all_IDs;
% end
