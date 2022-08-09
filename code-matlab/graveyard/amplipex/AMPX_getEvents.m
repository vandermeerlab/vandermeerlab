function [EventTimes] = APX_getEvents(data,condition)
% EventTimes will take in the file, number and setup and output the
% EventTimes structure containing the feeder fires, witheld pellets, nose poke entry leading to a
% fire, nose exits following a reward, nosepokes that did not lead to a
% reward
%
% 'data' is the channel file from AMPX_loadData.


%% Initialize Everything
addpath('C:\Program Files\MATLAB\R2013a\toolbox\signal\signal')
rmpath('C:\Program Files\MATLAB\R2013a\toolbox\chronux_2_10\chronux\spectral_analysis\continuous')
if strcmp(condition,'normal')== 1
    no_reward = 3; % this deals with which sites are rewarded under the normal and reversed condition
else
    no_reward = 5;
end
round = 1;  % %This sets the feeder id for each round

% if you want a graph for sanity checking purposes change graphs ==1
graphs = 1;
sub = 0; % the number of the subplot
%initialize all the Event ids and times
Events.fires.ids = [];             Events.fires.times = [];
Events.no_fires.ids= [];           Events.no_fires.times = [];
Events.nose_in_fires.ids = [];     Events.nose_in_fires.times = [];
Events.nose_out_fires.ids = [];    Events.nose_out_fires.times = [];
Events.nose_in_nofires.ids = [];   Events.nose_in_nofires.times = [];
Events.nose_out_nofires.ids = [];  Events.nose_out_nofires.times = [];

%% Initialize the channel and names
for kk = 2:5
    feeder = double(data.channels{kk});
    nose = double(data.channels{kk+4});
    withheld = double(data.channels{1});
    
    
    %% Identifies the feeder fires, nose-poke times and withheld fires
    [~,feeder_points] = findpeaks(feeder,'MINPEAKHEIGHT', 4000, 'MINPEAKDISTANCE', 6000);
    [~,withheld_points] = findpeaks(withheld,'MINPEAKHEIGHT', 4000, 'MINPEAKDISTANCE', 6000);
    [~,nose_points] = findpeaks(abs(nose),'MINPEAKHEIGHT', 4000, 'MINPEAKDISTANCE', 6000);
    %% Determines the feeder fires
    f_fires = nan(length(feeder_points),1);
    if kk ~= no_reward
        for jj = 1:length(feeder_points)
            f_fires(jj) = feeder_points(jj);
        end
    else
        continue
    end
    %% finds all the withheld trials.
    no_fires = nan(length(withheld_points),1);
    if kk ~= no_reward
        for jj = 1:length(withheld_points)
            no_fires(jj) = withheld_points(jj);
        end
    end
    
    %% Determines the nosepokes surrounding the feeder fires
    nose_enter_fires = NaN(length(feeder_points),1);
    nose_enter_fires = NaN(length(feeder_points),1);
    if kk ~= no_reward;
        for ii = 1:length(feeder_points)
            temp = feeder_points(ii) - nose_points;
            nosein_point = temp(6500>temp & temp>1000);
            if isempty(nosein_point) ==1
                disp('the delay was too short something is wrong')
                exception = temp < 1000;
                for jj = 1:length(exception)
                    if exception(jj) == 0 && exception(jj+1) == 1
                        nosein_point = temp(jj);
                        continue
                    end
                end
                nosein_loc = find(temp==nosein_point);
                nosein = nose_points(nosein_loc);
                noseout = nose_points(nosein_loc+1);
                nose_enter_fires(ii) = nosein;
                nose_exit_fires(ii) = noseout;
            else
                nosein_loc = find(temp==nosein_point);
                nosein = nose_points(nosein_loc);
                noseout = nose_points(nosein_loc+1);
                nose_enter_fires(ii) = nosein;
                nose_exit_fires(ii) = noseout;
            end
        end
    end
    %% Determines the nosepokes for withheld fires
    nose_enter_nofires = NaN(length(withheld_points),1);
    nose_exit_nofires = NaN(length(withheld_points),1);
    if kk ~= no_reward;
        for ii = 1:length(withheld_points)
            temp = withheld_points(ii) - nose_points;
            nosein_point = temp(6500>temp & temp>1000);
            if ~isempty(nosein_point)
                nosein_loc = find(temp==nosein_point);
                nosein = nose_points(nosein_loc);
                noseout = nose_points(nosein_loc+1);
                nose_enter_nofires(ii) = nosein;
                nose_exit_nofires(ii) = noseout;
            else
                disp('the delay was too short something is wrong')
                exception = temp < 1000;
                for jj = 1:length(exception)
                    if exception(jj) == 0 && exception(jj+1) == 1
                        nosein_point = temp(jj);
                        continue
                    end
                end
                nosein_loc = find(temp==nosein_point);
                nosein = nose_points(nosein_loc);
                noseout = nose_points(nosein_loc+1);
                nose_enter_nofires(ii) = nosein;
                nose_exit_nofires(ii) = noseout;
            end
        end
    end
    %% Plots all four analog channels with their nosepokes
    if graphs ==1
        sub = sub+1;
        subplot(4,1,sub)
        plot(1:5*10^6,nose(1:5*10^6),'b',1:5*10^6,feeder(1:5*10^6),'r')
    end
    
    Events.fires.times = [Events.fires.times , f_fires'];
    Events.fires.ids = [Events.fires.ids , (round*ones(1,length(f_fires)))];
    Events.no_fires.times = [Events.no_fires.times , no_fires'];
    Events.no_fires.ids = [Events.no_fires.ids , (round*ones(1,length(no_fires)))];
    Events.nose_in_nofires.times = [Events.nose_in_nofires.times , nose_enter_nofires'];
    Events.nose_in_nofires.ids = [Events.nose_in_nofires.ids , (round*ones(1,length(nose_enter_nofires)))];
    Events.nose_out_nofires.times = [Events.nose_out_nofires.times , nose_exit_nofires'];
    Events.nose_out_nofires.ids = [Events.nose_out_nofires.ids , (round*ones(1,length(nose_exit_nofires)))];
    Events.nose_in_fires.times = [Events.nose_in_fires.times , nose_enter_fires'];
    Events.nose_in_fires.ids = [Events.nose_in_fires.ids , (round*ones(1,length(nose_enter_fires)))];
    Events.nose_out_fires.times = [Events.nose_out_fires.times , nose_exit_fires];
    Events.nose_out_fires.ids = [Events.nose_out_fires.ids , (round*ones(1,length(nose_exit_fires)))];
    
    round = round+1;  % this is used for naming the feeder ids
end
save('EventData.mat','Events')
