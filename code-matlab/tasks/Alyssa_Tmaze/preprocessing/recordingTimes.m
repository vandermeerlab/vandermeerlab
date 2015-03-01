function rec_iv_pairs = recordingTimes
%RECORDINGTIMES Return the recording time intervals in start-stop pairs
%   For T-maze ExpKeys fields prerecord, task, postrecord
%
%   If more than 3 pairs exit, likely the last 3 pairs are the ones to be
%   kept, and the data needs to be looked into.
%
% A.Carey

evts = LoadEvents([]);


rec_start = getd(evts,'Starting Recording')';
rec_stop = getd(evts,'Stopping Recording')';

rec_iv_pairs = [rec_start rec_stop];

end

