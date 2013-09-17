function [EventTimes] = APX_getEvents(fname, rat_num,condition)
% EventTimes will take in the file, number and setup and output the
% EventTimes structure containing the feeder fires, witheld pellets, nose poke entry leading to a
% fire, nose exits following a reward, nosepokes that did not lead to a
% reward
%
% 'fname' is the name of the recording session as a string (ex:
% 'R041-2013-08-06'.  'rat_num' is the number of the rat (ex: 41).  The
% 'condition' is either 'normal' for the standard feeder setup or 'reversed'
% for the reversed feeder setup


%% Initialize Everything
addpath('C:\Program Files\MATLAB\R2013a\toolbox\signal\signal')
rmpath('C:\Program Files\MATLAB\R2013a\toolbox\chronux_2_10\chronux\spectral_analysis\continuous')
cd(['D:\DATA\R0' num2str(rat_num) '\' fname]);

for kk = 66:69
%% Initializes the names and calls
%     a = a+1;
    chan_call = ['channel_' num2str(kk)];
    chan_name = [fname '_ch' num2str(kk) '_2kHz.mat'];
    load(chan_name);
    feeder = chan;
    clear chan
    chan_name = [fname '_ch' num2str(kk+4) '_2kHz.mat'];
    load(chan_name);
    nose = chan;
    clear chan
    chan_name = [fname '_ch' num2str(65) '_2kHz.mat'];
    load(chan_name);
    withheld = chan;
    clear chan
    
%% Identifies the feeder fires, nose-poke times and withheld fires
    [~,feeder_points] = findpeaks(feeder,'MINPEAKHEIGHT', 4000, 'MINPEAKDISTANCE', 6000);
    [~,withheld_points] = findpeaks(withheld,'MINPEAKHEIGHT', 4000, 'MINPEAKDISTANCE', 6000);
    [~,nose_points] = findpeaks(abs(nose),'MINPEAKHEIGHT', 4000, 'MINPEAKDISTANCE', 6000);
%% Determines the feeder fires
    f_fires = nan(length(feeder_points),1);
    if kk ~= 67
        for jj = 1:length(feeder_points)
            f_fires(jj) = feeder_points(jj);
        end
    else
        continue
    end
%% finds all the withheld trials.
    no_fires = nan(length(withheld_points),1);
    if kk ~= 67
        for jj = 1:length(withheld_points)
            no_fires(jj) = withheld_points(jj);
        end
    else
        %         f_fires(jj) = [];
    end
    
%% Determines the nosepokes surrounding the feeder fires
    nf_fires = cell(length(feeder_points),1);
    for ii = 1:length(feeder_points)
        temp = feeder_points(ii) - nose_points;
        nosein_point = temp(7500>temp & temp>1000);
        nosein_loc = find(temp==nosein_point);
        nosein = nose_points(nosein_loc);
        noseout = nose_points(nosein_loc+1);
        nf_fires{ii} = [nosein, noseout];
    end
%% Determines the withheld fires 
    for ii = 1:length(withheld_points)
        temp = withheld_points(ii) - nose_points;
        nosein_point = temp(4000>temp & temp>1000);
        if ~isempty(nosein_point)
            nosein_loc = find(temp==nosein_point);
            nosein = nose_points(nosein_loc);
            noseout = nose_points(nosein_loc+1);
            nf_nofires{ii} = [nosein,noseout];
        else
            nf_nofires{ii} = [];
            continue
        end
    end
    
    
    
end
