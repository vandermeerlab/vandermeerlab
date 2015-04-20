function [ light_idx ] = light_on( position_x,position_y,Timestamps,event_file,event_id )

% Sort EventStrings and make an Event struct (from Youki)
fn = FindFile(event_file);
[EVTimeStamps, ~, ~, ~, EventStrings, ~] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

ev = unique(EventStrings);
for i = 1:length(ev);
    Events(i).type = ev(i);
    Events(i).index = find(strcmp(EventStrings,Events(i).type));
    Events(i).event_times = EVTimeStamps(Events(i).index);
    n_idx = nearest_idx3(Events(i).event_times,Timestamps);
    Events(i).position_x = position_x(n_idx);
    Events(i).position_y = position_y(n_idx);
    Events(i).n_idx = n_idx;
end

% Get times where LED is flashing (Event strings 10,11). 
% NOTE: each index in Events is the time that the LED turns on. 
% It STAYS ON for 3 timestamps after!!! 
% So, we make a new index list that contains three timesteps after.

light_idx = Events(event_id).n_idx;
temp_idx = [light_idx;light_idx+1;light_idx+2;light_idx+3];
light_idx = sort(temp_idx);
light_idx = light_idx(light_idx < length(position_x));
end