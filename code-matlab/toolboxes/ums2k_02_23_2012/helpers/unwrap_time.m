function  spiketimes  = unwrap_time( t, tr, dur, spacing)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/08/2010
%
% unwrap_time - convert from within trial time to absolute time
%
% Usage:
%           spiketimes  = unwrap_time( t, tr, dur, spacing)%
% Description:
%    Takes a list of spike times relative to the beginning of each trial
% and converts it to absolute time within the session.  A short period of
% time, given by "spacing", is padded between each trial.  It is assumed
% that trials occur in order of their number.
%
% Input: 
%   t     - [1 X N] list of spike event times relative to start of trial
%   tr    - [1 X N] trial in which event ocurred.
%   dur   - [1 X no. of trials] duration of each trial (s)
%   spacing  - amount of time to pad between each trial (s)
%
% Output:
%   spiketimes - [1 x N] list  of spike event times relative to 1st trial
%
    %determine absolute time of begninning of each trial
    cumdur = cumsum([0 dur(1:end-1)]) + spacing*[0:(length(dur)-1)];
    
    % add trial start time as an offset to each spike time
    spiketimes = t + cumdur(tr);

