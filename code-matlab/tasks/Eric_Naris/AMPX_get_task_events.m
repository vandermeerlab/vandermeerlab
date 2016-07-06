function [ evts ] = AMPX_get_task_events(cfg_in, data_AMPX)
%% AMPX_get_task_events
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
%%
cfg_def.win_r_in = -1; % pre reward period
cfg_def.win_r_out = 3; % post reward period
cfg_def.win_a_in = -5;  % time block before the start of the reward in
cfg_def.debug = 0;
cfg = ProcessConfig2(cfg_def, cfg_in);
%% find the peaks based on the cfg.trialdef.event_type
warning off
% Find Omissions:
[~,omit_peaks] = findpeaks(data_AMPX.channels{65}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);

% find the individual feeders
% feeder #1
[~,feeder1_fires] = findpeaks(data_AMPX.channels{66}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% feeder #2
[~,feeder2_fires] = findpeaks(data_AMPX.channels{67}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% feeder #3
[~,feeder3_fires] = findpeaks(data_AMPX.channels{68}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% feeder #4
[~,feeder4_fires] = findpeaks(data_AMPX.channels{69}, 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 1000);
% evts.task.omissions.times = omit_peaks;

warning on
%% create an ouput structure
fire1 = [feeder1_fires, ones(length(feeder1_fires),1)];
fire2 = [feeder2_fires, ones(length(feeder2_fires),1)*2];
fire3 = [feeder3_fires, ones(length(feeder3_fires),1)*3];
fire4 = [feeder4_fires, ones(length(feeder4_fires),1)*4];

omit = [omit_peaks, ones(length(omit_peaks),1)*5];

% get the photobeam breaks to determine the time at the feeders
% evts.task.photo1.times = feeder1_photo;
% evts.task.photo2.times = feeder2_photo;
% evts.task.photo3.times = feeder3_photo;
% evts.task.photo4.times = feeder4_photo;

%%
[~,feeder1_photo] = findpeaks(abs(data_AMPX.channels{70}), 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 500);
[~,feeder2_photo] = findpeaks(abs(data_AMPX.channels{71}), 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 500);
[~,feeder3_photo] = findpeaks(abs(data_AMPX.channels{72}), 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 500);
[~,feeder4_photo] = findpeaks(abs(data_AMPX.channels{73}), 'MINPEAKHEIGHT', 1500, 'MINPEAKDISTANCE', 500);

photo1 = [feeder1_photo, ones(length(feeder1_photo),1)];
photo2 = [feeder2_photo, ones(length(feeder2_photo),1)*2];
photo3 = [feeder3_photo, ones(length(feeder3_photo),1)*3];
photo4 = [feeder4_photo, ones(length(feeder4_photo),1)*4];

all_peaks_photo(:,1) = [photo1(:,1); photo2(:,1); photo3(:,1); photo4(:,1)];
all_peaks_photo(:,2) = [photo1(:,2); photo2(:,2); photo3(:,2); photo4(:,2)];

[~, sort_idx] = sort(all_peaks_photo(:,1));
sort_photo = all_peaks_photo(sort_idx,:);

%% hack way through photobeam peaks to find the first break at a feeder and
% the last one before a break at the next feeder
evt_num = 1;
next = 0;
for ii = 1:length(sort_photo)
    if ii == 1
        photo(evt_num,1) = sort_photo(ii,1);
        photo(evt_num,3) = sort_photo(ii,2);
        hold_ids = sort_photo(ii,2);
    end
    
    if next == 1
        photo(evt_num,1) = sort_photo(ii,1);
        photo(evt_num,3) = sort_photo(ii,2);
        hold_ids = sort_photo(ii,2);
        next = 0;
    elseif sort_photo(ii,2) ~= hold_ids;
        photo(evt_num,2) = sort_photo(ii-1,1);
        photo(evt_num,3) = sort_photo(ii-1,2);
        evt_num = evt_num+1;
        next = 1;
    end
end
% remove last one incase the session ended during a photobreak
photo = photo(1:end-1,:);
%% plot
if cfg.debug
    c_ord = linspecer(4);
    close all
    figure(10)
    maximize
    hold on
    for iphoto = 1:4
        plot(data_AMPX.tvec(1:4000*200), data_AMPX.channels{69+iphoto}(1:4000*200), 'color', c_ord(iphoto,:))
    end
    plot(photo(:,1)/2000, 4500, 'color', c_ord(1, :),'marker', '.' )
    % plot(feeder2_photo(1:3)/2000, 4500, 'color', c_ord(2, :),'marker', '*' )
    % plot(feeder3_photo(1:3)/2000, 4500, 'color', c_ord(3, :),'marker', '*' )
    % plot(feeder4_photo(1:3)/2000, 4500, 'color', c_ord(4, :),'marker', '*' )
    
    plot(photo(:,2)/2000, 4500, 'color', c_ord(1, :),'marker', '+' )
    % plot(feeder2_photo_out(1:3)/2000, 4500, 'color', c_ord(2, :),'marker', '+' )
    % plot(feeder3_photo_out(1:3)/2000, 4500, 'color', c_ord(3, :),'marker', '+' )
    % plot(feeder4_photo_out(1:3)/2000, 4500, 'color', c_ord(4, :),'marker', '+' )
    
    legend('70', '71', '72', '73')
    xlim([1 300])
end
%% collect the events
evts.task.feeder.idx = [fire1(:,1); fire2(:,1); fire3(:,1); fire4(:,1)];
evts.task.feeder.ids = [fire1(:,2); fire2(:,2); fire3(:,2); fire4(:,2)];

evts.task.omit.idx = omit(:,1);
evts.task.omit.ids = omit(:,2);
%% display stats
fprintf(['\nTask Events:\nFeeder#1 fires: ' num2str(length(feeder1_fires)) '\nFeeder#2 fires: ' num2str(length(feeder2_fires)) '\nFeeder#3 fires: ' num2str(length(feeder3_fires)) '\nFeeder#4 fires: ' num2str(length(feeder4_fires)) '\nOmissions: ' num2str(length(omit_peaks)) '\n'])

%% extract the evts as times periods around the fire/omit as "reward" periods

evts.feeder_iv.tstart = data_AMPX.tvec(evts.task.feeder.idx)+cfg.win_r_in;
evts.feeder_iv.tend = data_AMPX.tvec(evts.task.feeder.idx)+cfg.win_r_out;
evts.feeder_iv.type = 'iv'; evts.feeder_iv.firstTimestamp = 0; evts.feeder_iv.cfg.history.cfg = []; evts.feeder_iv.cfg.history.mfun = [];

evts.omit_iv.tstart = data_AMPX.tvec(evts.task.omit.idx)+cfg.win_r_in;
evts.omit_iv.tend = data_AMPX.tvec(evts.task.omit.idx)+cfg.win_r_out;
evts.omit_iv.type = 'iv'; evts.omit_iv.firstTimestamp = 0; evts.omit_iv.cfg.history.cfg = []; evts.omit_iv.cfg.history.mfun = [];

%% extract the 2 secons before the feeder
evts.approach_iv.tstart = data_AMPX.tvec(evts.task.feeder.idx)+cfg.win_a_in;
evts.approach_iv.tend = data_AMPX.tvec(evts.task.feeder.idx)+cfg.win_r_in;
evts.approach_iv.type = 'iv'; evts.approach_iv.firstTimestamp = 0; evts.approach_iv.cfg.history.cfg = []; evts.approach_iv.cfg.history.mfun = [];

evts.approach_omit_iv.tstart = data_AMPX.tvec(evts.task.omit.idx)+cfg.win_a_in;
evts.approach_omit_iv.tend = data_AMPX.tvec(evts.task.omit.idx)+cfg.win_r_in;
evts.approach_omit_iv.type = 'iv'; evts.approach_omit_iv.firstTimestamp = 0; evts.approach_omit_iv.cfg.history.cfg = []; evts.approach_omit.cfg.history.mfun = [];

%% get the IVs for the photobeam breaks.

evts.photo_iv.tstart = data_AMPX.tvec(photo(:,1));
evts.photo_iv.tend = data_AMPX.tvec(photo(:,2));
evts.photo_iv.type = 'iv'; evts.photo_iv.firstTimestamp = 0; evts.photo_iv.cfg.history.cfg = []; evts.photo_iv.cfg.history.mfun = [];

