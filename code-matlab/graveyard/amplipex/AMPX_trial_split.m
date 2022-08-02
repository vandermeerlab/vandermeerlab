function [data_trl, trials_def, times_cfg] = AMPX_trial_split(data_ft, times, WindowSize)
%AMPX_trial_split (data_ft, times, WindowSize)
%
% This will take the center of the gamma events and create an event with a
% start and end equal to the center time +/- the WindowSize. 
% INPUTS:
%
% data_ft: FT data structure
% times: [n x 2 array] these are the start and end times of the gamma
% events from the gamma detection program.  
%  WindowSize: value of the window lenght on either side of the center of
%  the gamma event.  0.05 is recommended for gamma events.  
%
% OUTPUTS:
%
% data_trl: [1 x 1 struct] FT strcutre with subfields
%   .trials {trial_num} (#chan,:)
%   .time 

trl_ctr = round(mean(times,2));

trl_leftwin = WindowSize; trl_left_s = trl_leftwin * data_ft.hdr.Fs;
trl_rightwin = WindowSize; trl_right_s = trl_rightwin * data_ft.hdr.Fs;

trl_start_s = 0; trl_len_s = repmat(trl_start_s,size(trl_ctr));

cfg = [];
cfg.trl = cat(2,trl_ctr-trl_left_s,trl_ctr+trl_right_s,trl_len_s);
times_cfg = cfg.trl;

data_trl = ft_redefinetrial(cfg,data_ft);
trials_def = 1:length(data_trl.trial);

end

