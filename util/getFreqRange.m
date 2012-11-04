function fid = getFreqRange(freqRange,F)
% function fid = getFreqRange(freqRange,F)
%
% returns indices in F (as returned by Chronux) to match input

[val,low_ind] = min(abs(F-freqRange(1)));
[val,hi_ind] = min(abs(F-freqRange(end)));
fid = low_ind:hi_ind;