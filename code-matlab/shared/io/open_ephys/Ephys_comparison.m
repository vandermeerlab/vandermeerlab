function [comp_data, cfg_out] = Ephys_comparison(Rat_id, sess_id1, sess_id2, cfg)
%% Ephys_comparison: a fast method for looking at two pre-processed acute 
%  recording sessions in order to see differences across conditions or 
%  recording depths
%
% Inputs: 
%  Rat_id {string}: rat number (eg: 'R055)
%  sess_id {String}: the session number aka the recording number (eg: "66")
%  cfg  [struct]: contains the output from the Ephys_preprocess as well as
%  any additional information such as the recording potions/ISI/intensity
%
%
%% Define the defaults
cd(['G:\Acute\' Rat_id '\all_data'])
fprintf(['Moving to ' Rat_id '\n'])
% find the data struct for data set 1
cfg.fname1 = [Rat_id '_' sess_id1];
data_1dir = [cfg.fname1 '*'];

data_1dir = FindFiles(data_1dir);
if isempty(data_1dir); error('Data1 file not found.  Be sure that the data has been preprocessed using Ephys_acute_preprocess'); end
for el = 1:length(data_1dir)
    if strcmp(data_1dir{el}(end-2:end), 'mat')
        cfg.fname1 = data_1dir{el};
    end
end

% find the data struct for data set 2
cfg.fname2 = [Rat_id '_' sess_id2];
data_2dir = [cfg.fname2 '*'];

data_2dir = FindFiles(data_2dir);
if isempty(data_2dir); error('Data2 file not found.  Be sure that the data has been preprocessed using Ephys_acute_preprocess'); end
for el = 1:length(data_2dir)
    if strcmp(data_2dir{el}(end-2:end), 'mat')
        cfg.fname2 = data_2dir{el};
    end
end
% load the data sets
data_1 = load(cfg.fname1);

data_2 = load(cfg.fname2);

%% plot the outputs of the two conditions
figHandles = get(0,'Children');
if isempty(figHandles) ==0 || sum(figHandles == 4000) > 0
    close(4000)
end
figure(4000)


