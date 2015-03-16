%% load behavioural data
load('C:\Users\mvdmlab\Dropbox\controlscripts\Jimmie\2015-02-17\R057-2015-02-17-light_behaviour.mat')
TrialInfo = TrialFnc(eventLog); % generate behavioural data

%% load AMPX maze data
fname = 'R057-2015-02-17-light.dat';
channels_to_load = 65:76; %non neural signals
cd('E:\Jimmie\R057\R057-2015-02-17');
data = AMPX_loadData(fname,channels_to_load);

%% Find time offset between matlab and AMPX **maybe switch to average instead of first trial**
threshold = 3000; %threshold for an event
refractory_period = 5; %period of time before passing threshold considered a new event
trackPB_1 = find(data.channels{10}> threshold); %find possible events
trackPB_1_new = trackPB_1(1);
for i = 2:length(trackPB_1)
    if data.tvec(trackPB_1(i))- data.tvec(trackPB_1(i-1)) > refractory_period %new event if far away enough in time
    trackPB_1_new(end+1) = trackPB_1(i);
    end
end

offsetT = data.tvec(trackPB_1_new(1)) - TrialInfo.trialT(1); %time offset between first trial initiation in MATLAB and AMPX

%% Load Spikes        
%cd('E:\Jimmie\R057\R057-2015-02-17');    
fc = FindFiles('*.t');
S = LoadSpikes(fc);

%% SDF
global tbin_centers
t = [0 800];
iC = 9;
binsize = 0.001; 
tbin_edges = t(1):binsize:t(2);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
 
spk_count = histc(Data(S{iC}),tbin_edges);
spk_count = spk_count(1:end-1);

binsize = 0.001; % in seconds
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.02./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
figure(20)
plot(tbin_centers,gau_sdf);
title('SDF'); ylabel('firing rate'); xlabel('time(s)');

%% Generate event plot for trials
global dataPoint %global so can interact with enscroll
for i = 1:length(TrialInfo.trialT) 
ZeroT_trials(i) = TrialInfo.trialT(i) + offsetT; %apply offset to MATLAB data
dataPoint(i) = find(tbin_centers > ZeroT_trials(i), 1, 'first'); %find times for all events
end
global event_number
global h1
event_number = 1;
figure(21);
plot(tbin_centers,gau_sdf);
set(gca,'XLim',[tbin_centers(dataPoint(event_number)-5000) tbin_centers(dataPoint(event_number)+10000)]);
hold on; h1 = plot(tbin_centers(dataPoint(event_number)),0:1:70,'color','black');
title('event plot by trial'); ylabel('firing rate'); xlabel('time(s)');

enscroll %enable scrolling through events

%% average response across trials
temp = -5000; %starting point in time for averaged-response
for ij = 1:length(dataPoint)
for ii = 1:15001 %collect data for -5 to +10 surrounding all events
all_trials(ii,ij) = gau_sdf(dataPoint(ij) + temp);
temp = temp + 1;
end
temp = -5000;
end

for jj = 1:15001 %average responses surorunding events
average_trials(jj,1) = mean(all_trials(jj,:));
end

figure(23);
subplot(3,2,1);
plot(-5:.001:10,average_trials);
hold on; plot(0,0:.05:5,'color','black');
title('trial-averaged response'); ylabel('firing rate'); ylim([0 8]);

%% rewarded vs unrewarded trials
temp = -5000;
rew_trial_num = 1;
unrew_trial_num = 1;
for ik = 1:length(TrialInfo.trialT)
    switch TrialInfo.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
        for jk = 1:15001
            rew_trials(jk,rew_trial_num) = gau_sdf(dataPoint(ik) + temp);
            temp = temp + 1;            
        end
        rew_trial_num = rew_trial_num + 1;
        case 0
        for jk = 1:15001
            unrew_trials(jk,unrew_trial_num) = gau_sdf(dataPoint(ik) + temp);
            temp = temp + 1;            
        end
        unrew_trial_num = unrew_trial_num + 1;
    end
    temp = -5000;
end

for jj = 1:15001 %generate averaged responses
avg_rew_trials(jj,1) = mean(rew_trials(jj,:));
avg_unrew_trials(jj,1) = mean(unrew_trials(jj,:));
end

figure(23);
subplot(3,2,3); plot(-5:.001:10,avg_rew_trials);
hold on; plot(0,0:.05:5,'color','black');
title('rewarded trial-averaged response'); ylabel('firing rate'); ylim([0 8]);
subplot(3,2,5); plot(-5:.001:10,avg_unrew_trials);
hold on; plot(0,0:.05:5,'color','black');
title('unrewarded trial-averaged response'); ylabel('firing rate'); xlabel('time(s)'); ylim([0 8]);

%% Generate event plot for nosepokes
global dataPoint_nosepokes %global so can interact with enscroll_nosepokes
nosepoke_counter = 1;
for i = 1:length(TrialInfo.nosepokeT)
    if TrialInfo.nosepokeT(i) ~= 0
ZeroT_nosepokes(i) = TrialInfo.nosepokeT(i) + offsetT; %apply offset to MATLAB data
dataPoint_nosepokes(nosepoke_counter) = find(tbin_centers > ZeroT_nosepokes(i), 1, 'first'); %find times for all events
nosepoke_counter = nosepoke_counter + 1;
    end
end
global event_number_nosepokes
global h2
event_number_nosepokes = 1;
figure(22);
plot(tbin_centers,gau_sdf);
set(gca,'XLim',[tbin_centers(dataPoint_nosepokes(event_number_nosepokes)-5000) tbin_centers(dataPoint_nosepokes(event_number_nosepokes)+10000)]);
hold on; h2 = plot(tbin_centers(dataPoint_nosepokes(event_number_nosepokes)),0:1:70,'color','black');
title('event plot by nosepoke'); ylabel('firing rate'); xlabel('time(s)');

enscroll_nosepokes %enable scrolling through events

%% average response across nosepokes
temp = -5000; %starting point in time for averaged-response
for ij = 1:length(dataPoint_nosepokes)
for ii = 1:15001 %collect data for -5 to +10 surrounding all nosepokes
all_nosepokes(ii,ij) = gau_sdf(dataPoint_nosepokes(ij) + temp);
temp = temp + 1;
end
temp = -5000;
end

for jj = 1:15001 %average responses surorunding nosepokes
average_nosepokes(jj,1) = mean(all_nosepokes(jj,:));
end

figure(23);
subplot(3,2,2);
plot(-5:.001:10,average_nosepokes);
hold on; plot(0,0:.05:5,'color','black');
title('nosepoke-averaged response'); ylim([0 8]);

%% rewarded vs unrewarded nosepokes
temp = -5000;
rew_num = 1;
unrew_num = 1;
for ik = 1:length(dataPoint_nosepokes)
    switch TrialInfo.rewarded(ik) %add to rewarded or unrewarded count depending on if trial was rewarded or not
        case 1
        for jk = 1:15001
            rew_nosepokes(jk,rew_num) = gau_sdf(dataPoint_nosepokes(ik) + temp);
            temp = temp + 1;            
        end
        rew_num = rew_num + 1;
        case 0
        for jk = 1:15001
            unrew_nosepokes(jk,unrew_num) = gau_sdf(dataPoint_nosepokes(ik) + temp);
            temp = temp + 1;            
        end
        unrew_num = unrew_num + 1;
    end
    temp = -5000;
end

for jj = 1:15001 %generate averaged responses
avg_rew_nosepokes(jj,1) = mean(rew_nosepokes(jj,:));
avg_unrew_nosepokes(jj,1) = mean(unrew_nosepokes(jj,:));
end

figure(23);
subplot(3,2,4); plot(-5:.001:10,avg_rew_nosepokes);
hold on; plot(0,0:.05:5,'color','black');
title('rewarded nosepoke-averaged response'); ylim([0 8]);
subplot(3,2,6); plot(-5:.001:10,avg_unrew_nosepokes);
hold on; plot(0,0:.05:5,'color','black');
title('unrewarded nosepoke-averaged response'); xlabel('time(s)'); ylim([0 8]);

%% save figures
saveas(figure(20),'E:\Jimmie\Jimmie\Analysis\R057\R057-2015-02-17-SDF.jpeg');
saveas(figure(21),'E:\Jimmie\Jimmie\Analysis\R057\R057-2015-02-17-trial.jpeg');
saveas(figure(22),'E:\Jimmie\Jimmie\Analysis\R057\R057-2015-02-17-nosepoke.jpeg');
saveas(figure(23),'E:\Jimmie\Jimmie\Analysis\R057\R057-2015-02-17-avg.jpeg');