function [events] = Ephys_event_check()
%% checks the cd for an all_channels.events file, and if it is found it will load it and give basic stats on the interval of the events and their duration.  

%% check for an events file
all_chan_loc = FindFiles('all_channels.events', 'StartingDirectory', '.', 'CheckSubdirs',0);
if isempty(all_chan_loc) == 1
    error('No all_channel.events file exists in cd')
    return
end
disp('File Found')
%% load the events file
[data, time, x] = load_open_ephys_data('all_channels.events');

%% check if there is anything in it
if isempty(data) || isempty(time) 
    warning('File contains no events')
    events = [];
    
else
%% print stats
fprintf('\n')
disp('Event Stats')
fprintf(['Inter Event interval: ', num2str(time(3)-time(1)) '\n'])
fprintf(['Event Duration:       ', num2str(time(2)-time(1)) '\n'])

%% give an output of a everything

events.time = time;
events.header = x;
events.channel_num = data;
end