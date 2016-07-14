function AMPX_Naris_fig_4_CSD(cfg_in)
%% AMPX_Naris_fig_4_kCSD: uses the csd output from AMPX_kCSD_example and
% AMPX_cycle_kCSD to produce figure 4 for Carmichael, Gmaz and van der
% Meer.  
%   
% 
% 
%          Inputs: 
%           - 
%           - 
%           - 
%          Outputs: 
%           - 
%           - 
%           -  
% 
% EC - 2016-07-09

%% Collect inputs/defaults
cfg_def.session_name = 'R061_2014_09_26';
cfg_def.example = 29;

cfg = ProcessConfig(cfg_def, cfg_in)


%% generate a CSD for the same event used in figure 1c
% load('C:\temp\Naris_all_data_pre.mat');
example_out = AMPX_kCSD_example(cfg, all_data_pre)


%%