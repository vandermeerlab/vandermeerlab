function Ephys_chan_plot(data, cfg)
% Ephys_chan_plot(data, cfg): plots all the channels that have been loaded
% using Ephys_load_data.  Labels are applied based on the channel names.  x
% range can be specified in the cfg.
%
%inputs:
%   data [struct]: output from Ephys_load_data
%              - data
%              - tvec
%              - labels
%              - hdr (Fs, ...)
%   cfg [struct]: input parameters file

%% define defaults
if isfield(cfg, 'chan_to_view') ==0
    cfg.chan_to_view = 1:length(data.Channels); 
    disp(['Channels set to default: ' num2str(1:length(data.Channels))])
end
if isfield(cfg, 'offset') ==0
    cfg.offset = 500; disp(['Offset set to default: ' num2str(cfg.offset)])
end



%% Plot the channels of interest
color = linspecer(length(cfg.chan_to_view));
figHandles = get(0,'Children');
if isempty(figHandles) ==0 || sum(figHandles == 10) > 0
    close(10)
end
figure(10)
hold on
loop= 1;
for iChan = 1:length(cfg.chan_to_view)
    plot(((1:length(data.Channels{iChan}.data)).*data.Channels{iChan, 1}.info.header.sampleRate)/60, data.Channels{iChan}.data+cfg.offset*(-loop+1), 'Color', color(loop,:))
    leg_names{loop} =  data.Channels{iChan}.info.header.channel;
    if isfield(cfg, 'events_pos') ==1
            plot(cfg.events_pos*data.Channels{iChan, 1}.info.header.sampleRate/60, cfg.events_pow, 'rx')
    end
    loop = loop +1;
end
legend(leg_names)