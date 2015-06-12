function [data, cfg] = Ephys_load_session(data_dir, fname, cfg)
% Ephys_load_data: loads the data channels for the continuous recordings
% for the data specified in the fname folder.  Assumes only one recording
% exists in each folder.  (needs a way to adapt to multiple recordings in a
% single recording session)


%% 
tic
title_dir = fname; 
cfg.title_dir = title_dir;
cd([data_dir])
data_dir = [fname '*'];
data_dir = FindFiles(data_dir);
data_dir = data_dir{1,1};
disp(['Folder to load: ' fname])
cd(data_dir);
listings = dir;
if isempty(FindFiles('115*'))==0 && isempty(FindFiles('110*')) && isempty(FindFiles('100*'))
    prefix = '115_CH';
elseif isempty(FindFiles('110*')) ==0  && isempty(FindFiles('100*'))
    prefix = '110_CH';
elseif isempty(FindFiles('110*')) ==1  && isempty(FindFiles('100*'))==0
    prefix = '100_CH';
end

%% load the data
loop_num = 0;
for iChan = cfg.chan_to_view
    loop_num = loop_num+1;
    [data.Channels{loop_num, 1}.data,data.Channels{loop_num, 1}.tvec,data.Channels{loop_num, 1}.info] = load_open_ephys_data([prefix num2str(iChan) '.continuous']);
    data.Channels{loop_num, 1}.data = decimate(data.Channels{loop_num, 1}.data, cfg.decimate_factor);
    data.Channels{loop_num, 1}.tvec = decimate(data.Channels{loop_num, 1}.tvec, cfg.decimate_factor);
    data.Channels{loop_num, 1}.info.header.sampleRate = data.Channels{loop_num, 1}.info.header.sampleRate/cfg.decimate_factor;
    data.Channels{loop_num, 1}.label = iChan;
end
disp(['Channels loaded: ' num2str(cfg.chan_to_view)])
data.labels = cfg.chan_to_view;
toc