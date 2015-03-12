%% Open Ephys Data Analysis Sandbox

data.info = get_session_info(pwd);
% info.settings = xmlread('settings.xml');
num_chan = length(data.info.electrodes);
num_chan = 32

%% load all the data files.
data.Channels = cell(num_chan, 1);

for iChan = num_chan:-1:1
    [data.Channels{iChan}.data,data.Channels{iChan}.tvec,data.Channels{iChan}.info] = load_open_ephys_data(['100_CH' num2str(iChan) '.continuous']);
end

%% Try ploting some channels

chan_to_plot = 1:ceil(num_chan/8):num_chan;
% chan_to_plot = [9 12 17];
color = linspecer(length(chan_to_plot));
figure
hold on
% range = [data.Channels{iChan}.tvec(1) data.Channels{iChan}.tvec(30000/2)];
offset = 200;
round= 1;
for iChan = 1:length(chan_to_plot)
    plot(data.Channels{iChan}.tvec, data.Channels{iChan}.data+offset*round, 'Color', color(round,:))
    round = round +1;
end

% xlim(range)
% xlim([264 266])

%%
figure 
plot(data.Channels{10,1}.tvec, data.Channels{10,1}.data+200, 'b')
hold on
plot(data.Channels{9,1}.tvec, data.Channels{9,1}.data, 'r')
