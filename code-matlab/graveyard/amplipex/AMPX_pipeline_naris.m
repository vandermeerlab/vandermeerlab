function [data_out] = AMPX_pipeline_sfn(file_name, varargin)
%% AMPX_pipeline performs: powerspectral densities, phase offset and PCA on a given data session.
%                          This uses:
%                          AMPX_loadData
%                          AMPX_RemoveArtifacts
%                          AMPX_filter
%                          AMPX_trial_split
%                          AMPX_pow_distrib
%                          AMPX_pow_distrib_plot_naris
%                          AMPX_phase_distrib
%                          AMPX_phase_distrib_plot
%
%
%INPUTS:
%   - file_name: 'R0XX-yyyy-mm-dd' If this is a pre or post session add
%   "-post" or "-pre"
% Optional:
%   - 'gamma_type' : either 'high', 'low' or if left empty will do both
%   - 'figures' : 'on' or 'off' [default: 'off']
%   - 'savefig' : 'on' or 'off' [default: 'off']
%   - 'gamma_chan' : channel for gamma detection.  This is normally done
%   using the channel with the highest power.
%OUTPUTS:
%   - 'data_out' : a structure containing the power/phase/PCA results for
%   both the high and low gamma events as well as some of the parameters
%   used.
%
% EC 09/2014
%% Set some variables
save_fig = 'on';
figures = 'on';
% display = 'off'; % used for the gamma detection viewer
peri_gamma_evt_window = .25;
gamma_evt_window = 0.025;
% theta_evt_window = 0.5;
% delta_evt_window = 1.5;
decimate_factor = 10;
% event_ID = 1;
extract_varargin

%% determine the session type
file_name = strrep(file_name, '_', '-');
if isempty(str2num(file_name(end)))==1
    session_name = file_name(1:15);
else
    session_name = file_name;
end

cd(['D:\DATA\' session_name(1:4) '\' session_name]);
current_cd = cd;

if strcmp(file_name(end-3:end), 'post')
    session_type = 'post';
elseif strcmp(file_name(end-2:end), 'pre')
    session_type = 'pre';
else
    session_type = 'task';
end
%% load the data
if exist([session_type '_data_preprocess_filt.mat'],'file')==0
    if exist([ session_type '_data_preprocess.mat'],'file')~=0
        disp('Data file found')
        load([session_type '_data_preprocess.mat'])
        disp('Data file loaded')
    else
        disp('No data file found..')
        disp('Loading Data');
        % channels_to_load = [33,57,1,25,36,58,2,26,38,62,6,30,37,61,5,29]; % specify channels in the APX data file to load
        if strcmp(session_type, 'task') == 1;
            channels_to_load = 1:74;
        else
            channels_to_load = 1:64;
        end
        data = AMPX_loadData([file_name '.dat'],channels_to_load,decimate_factor); % note decimation factor of 20
        disp('Data loading complete');
        save([session_type '_data_preprocess.mat'], 'data', '-v7.3')
    end
    session_name = [data.hdr.Filename(5:8) '_' data.hdr.Filename(10:13) '_' data.hdr.Filename(15:16) '_' data.hdr.Filename(18:19)];
    
    %% DC remove (only if this has not previously been done.
    if isfield(data,'data_ft') ==0
        for iChan = size(data.channels,2):-1:1
            data_temp{1,iChan} = data.channels{1,iChan} - mean(data.channels{1,iChan});
        end
        data.channels = data_temp;
        data.dc_remove = 'yes';
        clear data_temp;

        %% Bandpass the signal
        [data.data_ft] = AMPX_filter(data, 1, 500, 10);
        
        save([session_type '_data_preprocess_filt.mat'], 'data', '-v7.3')
    end
else
    load([session_type '_data_preprocess_filt.mat'])
    session_name = [data.hdr.Filename(5:8) '_' data.hdr.Filename(10:13) '_' data.hdr.Filename(15:16) '_' data.hdr.Filename(18:19)];
end
%% pull the data type out of the data struct

data_ft = data.data_ft;
data = rmfield(data, 'data_ft');

%% theta and Spindle detection.

% [theta_iv, spindles_iv] = AMPX_theta_spindle_detect(data);


%% Delta detection

% [delta_iv] = AMPX_delta_detect(data);


%%  Load the all events file
if exist('D:\DATA\Si_probe_analysis_data\all_events.mat', 'file')
    load('D:\DATA\Si_probe_analysis_data\all_events.mat')
else %create a local temporary file
    fprintf('No all events file found\n')
    all.ev = [];
    fprintf('temporary all.ev file created in the CD\n')
end
%% Reformat the data to fit the task events.
if strcmp(session_type, 'task')==1
    load('events.mat')
    event_ID = 1:5;
    [all_events] = AMPX_format_task_events(data, evt);
    disp('Data has been recompiled based on task event segments')
end

%% Delta Detection
% mkdir([date '\' session_type '\Delta_detect_figures'])
% 
% all.ev.(session_name).(session_type).delta = ([delta_iv.tstart, delta_iv.tend])*data_ft.hdr.Fs;
% if isempty(delta_iv.tstart) ==1
%     data_out.delta.power.power_distrib = [];
% else
%     [data_trl_delta, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).delta, delta_evt_window);
%     [data_out.delta.power] = AMPX_pow_distrib(data_trl_delta, 'min_freq', 2, 'max_freq', 4);
% %         [data_out.delta.phase] = AMPX_phase_distrib(data_trl_delta, 'min_freq', 2, 'max_freq', 4);
%     if strcmp(figures, 'on')==1
%         mkdir([date '\' session_type '\delta_pow_figures']);  mkdir([date '\' session_type '\delta_phase_figures']);
%         AMPX_pow_distrib_plot_naris_avg(data_out.delta.power, data_trl_delta, data, 'save_fig', save_fig)
% %                 AMPX_phase_distrib_plot(data_out.theta.phase, data_trl_theta, data, 'save_fig', save_fig)
%     end
% end

%% gamma detection
mkdir([date '\' session_type '\Gamma_detect_figures'])

[low_gamma_iv, high_gamma_iv, random_lg_iv, random_hg_iv, ~, ~] = AMPX_Gamma_detect_nsb2014 (data, 'display', 'off');

all.ev.(session_name).(session_type).high_gamma = ([high_gamma_iv.tstart, high_gamma_iv.tend])*data_ft.hdr.Fs;
all.ev.(session_name).(session_type).low_gamma = ([low_gamma_iv.tstart, low_gamma_iv.tend])*data_ft.hdr.Fs;
% randoms
all.ev.(session_name).(session_type).random_high_gamma = ([random_hg_iv.tstart, random_hg_iv.tend])*data_ft.hdr.Fs;
all.ev.(session_name).(session_type).random_low_gamma = ([random_lg_iv.tstart, random_lg_iv.tend])*data_ft.hdr.Fs;
close all
disp('detection Complete')


%% section the events and filter for each bands of interest
% high gamma
if isempty(high_gamma_iv.tstart) ==1
    data_out.high_gamma.power.power_distrib = [];
else
    [data_trl_high, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).high_gamma, gamma_evt_window);
    [data_trl_high_peri, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).high_gamma, peri_gamma_evt_window);
    [data_out.high_gamma.power] = AMPX_pow_distrib(data_trl_high, 'min_freq', 70, 'max_freq', 85);
    [data_out.high_gamma.phase] = AMPX_phase_distrib(data_trl_high, 'min_freq', 70, 'max_freq', 85);
    
    if strcmp(figures, 'on')==1
        mkdir([date '\' session_type '\high_gamma_pow_figures']);  mkdir([date '\' session_type  '\high_gamma_phase_figures']);
%         AMPX_pow_distrib_plot_naris(data_out.high_gamma.power, data_trl_high_peri, data, 'save_fig', save_fig, 'fig_type', 'high')
        AMPX_pow_distrib_plot_naris_avg(data_out.high_gamma.power, data_trl_high_peri, data, 'save_fig', save_fig, 'fig_type', 'high')
%         AMPX_phase_distrib_plot(data_out.high_gamma.phase, data_trl_high, data, 'save_fig', save_fig, 'fig_type', 'high')
    end
end
%% Low Gamma
if isempty(low_gamma_iv.tstart) ==1
    data_out.low_gamma.power.power_distrib = [];
else
    [data_trl_low, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).low_gamma, gamma_evt_window);
    mkdir([date '\' session_type '\low_gamma_pow_figures']); mkdir([date '\' session_type '\low_gamma_phase_figures']);
    [data_out.low_gamma.power] = AMPX_pow_distrib(data_trl_low, 'min_freq', 40, 'max_freq', 55);
      tic; [data_out.low_gamma.phase] = AMPX_phase_distrib(data_trl_low, 'min_freq', 40, 'max_freq', 55); toc;
    
    if strcmp(figures, 'on')==1
        [data_trl_low_peri, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).low_gamma, peri_gamma_evt_window);
        AMPX_pow_distrib_plot_naris_avg(data_out.low_gamma.power, data_trl_low_peri, data, 'save_fig', save_fig,'fig_type', 'low')
%         AMPX_pow_distrib_plot_naris(data_out.low_gamma.power, data_trl_low_peri, data, 'save_fig', save_fig,'fig_type', 'low')
%         AMPX_phase_distrib_plot(data_out.low_gamma.phase, data_trl_low, data, 'save_fig', save_fig,'fig_type', 'low')
    end
end

%% random "low" intervals
if isempty(random_lg_iv.tstart) ==1
    data_out.random_low_gamma.power.power_distrib = [];
else
    [data_trl_low, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).random_low_gamma, gamma_evt_window);
    mkdir([date '\' session_type '\Ran_low_gamma_pow_figures']); mkdir([date '\' session_type '\Ran_low_gamma_phase_figures']);
    [data_out.random_low_gamma.power] = AMPX_pow_distrib(data_trl_low, 'min_freq', 40, 'max_freq', 55);
     tic; [data_out.random_low_gamma.phase] = AMPX_phase_distrib(data_trl_low, 'min_freq', 40, 'max_freq', 55); toc;
    
    if strcmp(figures, 'on')==1
        AMPX_pow_distrib_plot_naris_avg(data_out.random_low_gamma.power, data_trl_low, data, 'save_fig', save_fig,'fig_type', 'Ran_low')
%         AMPX_pow_distrib_plot_naris(data_out.random_low_gamma.power, data_trl_low, data, 'save_fig', save_fig,'fig_type', 'Ran_low')
%         AMPX_phase_distrib_plot(data_out.random_low_gamma.phase, data_trl_low, data, 'save_fig', save_fig, 'fig_type', 'Ran_low')
    end
end
%% random "high" intervals
if isempty(random_hg_iv.tstart) ==1
    data_out.random_high_gamma.power.power_distrib = [];
else
    [data_trl_high_ran, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).random_high_gamma, gamma_evt_window);
    mkdir([date '\' session_type '\Ran_high_gamma_pow_figures']); mkdir([date '\' session_type '\Ran_high_gamma_phase_figures']);
    [data_out.random_high_gamma.power] = AMPX_pow_distrib(data_trl_high_ran,  'min_freq', 70, 'max_freq', 85);
     tic; [data_out.random_high_gamma.phase] = AMPX_phase_distrib(data_trl_high_ran,  'min_freq', 70, 'max_freq', 85); toc;
    
    if strcmp(figures, 'on')==1
%         AMPX_pow_distrib_plot_naris(data_out.random_high_gamma.power, data_trl_high_ran, data, 'save_fig', save_fig,'fig_type', 'Ran_high')
        AMPX_pow_distrib_plot_naris_avg(data_out.random_high_gamma.power, data_trl_high_ran, data, 'save_fig', save_fig,'fig_type', 'Ran_high')
%         AMPX_phase_distrib_plot(data_out.random_high_gamma.phase, data_trl_high_ran, data, 'save_fig', save_fig,'fig_type', 'Ran_high')
    end
end
% %% Theta detection
% mkdir([date '\' session_type '\Theta_detect_figures'])
% 
% all.ev.(session_name).(session_type).theta = ([theta_iv.tstart, theta_iv.tend])*data_ft.hdr.Fs;
% if isempty(theta_iv.tstart) ==1
%     data_out.theta.power.power_distrib = [];
% else
%     [data_trl_theta, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).theta, theta_evt_window);
%     [data_out.theta.power] = AMPX_pow_distrib(data_trl_theta, 'min_freq', 6, 'max_freq', 8);
% %     [data_out.theta.phase] = AMPX_phase_distrib(data_trl_theta, 'min_freq', 7, 'max_freq', 10);
%     if strcmp(figures, 'on')==1
%         mkdir([date '\' session_type '\theta_pow_figures']);  mkdir([date '\' session_type '\theta_phase_figures']);
%         AMPX_pow_distrib_plot_naris_avg(data_out.theta.power, data_trl_theta, data, 'save_fig', save_fig)
% %         AMPX_phase_distrib_plot(data_out.theta.phase, data_trl_theta, data, 'save_fig', save_fig)
%     end
% end
%% Spindles
% mkdir([date '\' session_type '\spindles_detect_figures'])
% all.ev.(session_name).(session_type).spindles = ([spindles_iv.tstart, spindles_iv.tend])*data_ft.hdr.Fs;
% if isempty(spindles_iv.tstart) ==1
%     data_out.spindles.power.power_distrib = [];
% else
%     [data_trl_spindles, ~, ~] = AMPX_trial_split(data_ft, all.ev.(session_name).(session_type).spindles, theta_evt_window);
%     [data_out.spindles.power] = AMPX_pow_distrib(data_trl_spindles, 'min_freq', 8.5, 'max_freq', 11);
% %         [data_out.spindles.phase] = AMPX_phase_distrib(data_trl_spindles, 'min_freq', 8.5, 'max_freq', 10.5);
%     if strcmp(figures, 'on')==1
%         mkdir([date '\' session_type '\spindles_pow_figures']);  mkdir([date '\' session_type '\spindles_phase_figures']);
%         AMPX_pow_distrib_plot_naris_avg(data_out.spindles.power, data_trl_spindles, data, 'save_fig', save_fig)
% %                 AMPX_phase_distrib_plot(data_out.spindles.phase, data_trl_spindles, data, 'save_fig', save_fig)
%     end
% end
% %% Look at the difference in the plots for task related events.
% if strcmp(session_type, 'task')
%     events_out = AMPX_feeder_event_intersect(all_events, all.ev.(session_name).(session_type), data, 'pre_in', 2, 'post_win', 4);
%     data_out.events_out = events_out;
% end
%% Run PCA on all the bands of interest
% if isempty(high_gamma_iv.tstart) ~=1
%     [data_out.high_gamma.PCA] = AMPX_PCA(data_out.high_gamma.power, data);
% end
% if isempty(low_gamma_iv.tstart) ~=1
%     [data_out.low_gamma.PCA] = AMPX_PCA(data_out.low_gamma.power, data);
% end
% if isempty(theta_iv.tstart) ~=1
%     [data_out.theta.PCA] = AMPX_PCA(data_out.theta.power, data);
% end
% if isempty(spindles_iv.tstart) ~=1
%     [data_out.spindles.PCA] = AMPX_PCA(data_out.spindles.power, data);
% end
% if isempty(delta_iv.tstart) ~=1
%     [data_out.delta.PCA] = AMPX_PCA(data_out.delta.power, data);
% end
%% add in a cfg to the data_out
data_out.cfg.fname = session_name;
data_out.cfg.type = session_type;
data_out.cfg.version = '1.0';
%% Save the Data_out
 save([cd '\' date '\' session_type '_data_out2.mat'], 'data_out','-v7.3'); save([session_type '_data_out'], 'data_out','-v7.3');
% save 'D:\DATA\Si_probe_analysis_data\all_events.mat' all
fprintf(['\n' file_name ' Complete'])