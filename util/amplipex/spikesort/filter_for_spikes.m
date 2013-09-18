function out = filter_for_spikes(in,varargin)
% function out = filter_for_spikes(in,varargin)
%
%

Fs = 20000;
extract_varargin;

Wp = [ 700 7000] * 2 / Fs; % pass band for filtering
Ws = [ 500 9000] * 2 / Fs; % transition zone
[N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
[B,A] = butter(N,Wn); % builds filter
out = filtfilt( B, A, double(in) ); % runs filter
