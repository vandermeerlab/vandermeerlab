function ts_out = ts(varargin)
% function ts_out = ts(varargin)
%
% constructor for ts (timestamp data) struct
%
% MvdM 2014-06-17

ts_out.t = {};
ts_out.label = {};


% housekeeping
ts_out.cfg.history.mfun{1} = mfilename;
ts_out.cfg.history.cfg{1} = [];