function [Labels] = AMPX_ChannelLabels(data)

addpath('C:\Users\mvdmlab\Dropbox\Matlab\Analysis\Channel_labels')

%%load the Exp_keys to find out what the conditions for that session were
% if the session was a reversal or not.  Which HS was in which position.
% In which orientation the HS and preamps were connected.
fname = strrep(data.hdr.Filename, 'D://','');
fname = strrep(fname, '.dat','');
load(['fname','.m']);


