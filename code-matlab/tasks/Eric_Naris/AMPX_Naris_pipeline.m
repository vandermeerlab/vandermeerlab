function [data_out] = AMPX_Naris_pipeline(fname, varargin)
%% AMPX_Naris_pipeline performs: powerspectral densities, phase offset and PCA on a given data session.
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
% EC 04/2016


%% Preamble

% CD to the data directory.  Should contain all the .dat, .ini, and .meta
cd('D:\DATA\R054\R054-2014-10-10')
mkpre % gets the file name for the pre record.  
%%  determine session type
fname = strrep(file_name, '_', '-');
% is this sessions a pre, task or post record?
if strcmp(fname(end-3:end), 'post')
    session_type = 'post';
elseif strcmp(fname(end-2:end), 'pre')
    session_type = 'pre';
else
    session_type = 'task';
end

if isempty(str2num(fname(end)))==1
    fname = fname(1:15);
end

%% preprocess/load the data

[data, data_ft] = AMPX_Naris_preprocess([],file_name,session_type);

LoadExpKeys

%% get the events (using Catanese et al. method)

[evts, detect_csc] = AMPX_Julien_DetectEvents(data, ExpKeys); 

%% get an equal number of non-overlapping 
ran_evts = MatchGammaEvents([], detect_csc, evts.low, evts.high);
% evts = AMPX_get_random_evts(data, evts);
%% Organize the data to fit the probe lay out in the 8x8 format.
% Data about the probemapping 

data_remap_AMPX = AMPX_remap_sites(data, ExpKeys)

data_remap_FT = AMPX_remap_sites(data_ft, ExpKeys)

%% get the events (using Catanese et al. method)





