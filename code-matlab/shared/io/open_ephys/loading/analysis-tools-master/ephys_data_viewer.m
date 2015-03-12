function ephys_data_viewer(data_dir, chan_to_view, varargin)
%% ephys_data_viewer will import the continuous data files from the
% directory specified in 'data_dir'.  The data will then be decimated along
% with the tvec.  All these data sets will then be ploted.  The optional
% "peak_jumper" will create a subplot with the event centered and zoomed
% based on the "window" value.  Pressing any key will move through the
% events. Pressing 'ctrl + c' will exit the peak_jumper cycle.
%
%
%INPUTS:
%    - 'data_dir' {str}: path to file (eg:'G:\R047\2014-11-06_20-02-52')
%    - 'chan_to_view [1xN array]: list the channels to view (eg: [9 11])
%    - to use the other inputs type the name in '' then a comma, then the
%    value in the same format as it is listed in the Gather Variables
%    section.  Eg: if you want to change the decimation factor type it as:
%    ephys_data_viewer('G:\Acute\R047\HC_Stim_2014-12-04_14-15-33_r3', [9 11], 'decimate_factor', 3).
%    This can be stacked by putting commas between them Eg:
%    ephys_data_viewer('G:\Acute\R047\HC_Stim_2014-12-04_14-15-33_r3', [9 11], 'decimate_factor', 3, 'peak_jumper', 'on');

%% Gather the variables

decimate_factor = 10; % default is 3;
offset  = 500;       %used to separate the plots (default is 500)
peak_jumper = 'off'; % turn this 'on' to jump between events by hitting any key.
ISI = 10;            % Inter-Stimulus-interval in seconds, used for the peak jumper.
window = .5;        % display window in seconds, used for the peak jumper.
extract_varargin

cd(data_dir);


%% load the data, then decimate it
for iChan = chan_to_view
    [data.Channels{iChan, 1}.data,data.Channels{iChan, 1}.tvec,data.Channels{iChan, 1}.info] = load_open_ephys_data(['100_CH' num2str(iChan) '.continuous']);
    data.Channels{iChan,1}.data = data.Channels{iChan,1}.data(1:end-(round(.25*length(data.Channels{iChan,1}.data))));
    data.Channels{iChan,1}.tvec = data.Channels{iChan,1}.tvec(1:end-(round(.25*length(data.Channels{iChan,1}.tvec))));
    data.Channels{iChan, 1}.data = decimate(data.Channels{iChan, 1}.data, decimate_factor);
    data.Channels{iChan, 1}.tvec = decimate(data.Channels{iChan, 1}.tvec, decimate_factor);
    data.Channels{iChan, 1}.info.header.sampleRate = data.Channels{iChan, 1}.info.header.sampleRate/decimate_factor;
end

disp('Loading and Decimcation Complete')
%% Plot the channels of interest
color = linspecer(length(chan_to_view));
figure(10)
if strcmp(peak_jumper, 'on'); subplot(2,1,1) ; end
hold on
loop = 1;
for iChan = chan_to_view
    plot(1:length(data.Channels{iChan}.data), data.Channels{iChan}.data+offset*(-loop+1), 'Color', color(loop,:))
    leg_names{loop} =  data.Channels{iChan}.info.header.channel;
    loop = loop +1;
end
legend(leg_names)
%% subplot for visualizing the events one by one.
if strcmp(peak_jumper, 'on')
    color = linspecer(length(chan_to_view));
    subplot(2,1,2)
    hold on
    loop= 1;
    for iChan = chan_to_view
        %     plot(data.Channels{iChan}.tvec, data.Channels{iChan}.data+offset*(-loop+1), 'Color', color(loop,:))
        plot(1:length(data.Channels{iChan,1}.data), data.Channels{iChan}.data+offset*(-loop+1), 'Color', color(loop,:))
        leg_names{loop} =  data.Channels{iChan}.info.header.channel;
        loop = loop +1;
    end
    legend(leg_names)
    
    %% Peak jumper
    [~, peak_loc] = findpeaks(abs(data.Channels{chan_to_view(1)}.data), 'MinPeakHeight', 1000, 'MinPeakDistance', floor((data.Channels{iChan}.info.header.sampleRate)*(ISI-1)));
    for ii = 1:length(peak_loc)
        peak_x_val = data.Channels{chan_to_view(1)}.tvec(peak_loc(ii)) ;
        %         peak_vals(ii) = peak_x_val;
        xlim([peak_x_val-window peak_x_val+window]); % sets the x axis to +/- the window value around each peak.
        ylim([-loop*500 500])
        pause;
    end
end