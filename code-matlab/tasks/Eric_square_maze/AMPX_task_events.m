function [ evts ] = AMPX_task_events(data_filtered)
%% AMPX_task_events
% AMPX_task_events will use the analog inputs to extract events related to
% "nosepokes', 'ommissions', 'rewards'
%
%Site IDs: feeder 1 = 1; 
%          feeder 2 = 2; 
%          feeder 3 = 3; 
%          feeder 4 = 4; 
%          Omissions = 5; 
%
%OUTPUT: 
%  -evt [struct] contains:
%        - evt.task.times which contains all the event times in one array
%        - evt.task.times contains the corresponding ID numbers (see above)
%
%% find the peaks based on the cfg.trialdef.event_type
warning off
% Find Omissions:
[~,omit_peaks] = findpeaks(data_filtered.channels{65}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);

% find the individual feeders
% feeder #1
[~,feeder1_fires] = findpeaks(data_filtered.channels{66}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
[~,feeder1_photo] = findpeaks(data_filtered.channels{70}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% feeder #2
[~,feeder2_fires] = findpeaks(data_filtered.channels{67}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
[~,feeder2_photo] = findpeaks(data_filtered.channels{71}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% feeder #3
[~,feeder3_fires] = findpeaks(data_filtered.channels{68}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
[~,feeder3_photo] = findpeaks(data_filtered.channels{72}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% feeder #4
[~,feeder4_fires] = findpeaks(data_filtered.channels{69}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
[~,feeder4_photo] = findpeaks(data_filtered.channels{73}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);

warning on
%% create an ouput structure
% evts.task.omissions.times = omit_peaks;
fire1 = [feeder1_fires, ones(length(feeder1_fires),1)];
% evts.task.times = feeder1_photo;
fire2 = [feeder2_fires, ones(length(feeder2_fires),1)*2];
% evts.task.photo2.times = feeder2_photo;
fire3 = [feeder3_fires, ones(length(feeder3_fires),1)*3];
% evts.task.photo3.times = feeder3_photo;
fire4 = [feeder4_fires, ones(length(feeder4_fires),1)*4];

omit = [omit_peaks, ones(length(omit_peaks),1)*5];

% evts.task.photo4.times = feeder4_photo;
evts.task.times = [fire1(:,1); fire2(:,1); fire3(:,1); fire4(:,1); omit(:,1)];
evts.task.ids = [fire1(:,2); fire2(:,2); fire3(:,2); fire4(:,2); omit(:,2)];

%% display stats
fprintf(['\nTask Events:\nFeeder#1 fires: ' num2str(length(feeder1_fires)) '\nFeeder#2 fires: ' num2str(length(feeder2_fires)) '\nFeeder#3 fires: ' num2str(length(feeder3_fires)) '\nFeeder#4 fires: ' num2str(length(feeder4_fires)) '\nOmissions: ' num2str(length(omit_peaks)) '\n'])

end

