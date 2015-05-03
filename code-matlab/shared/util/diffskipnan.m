function out = diffskipnan(in)
% function out = diffskipnan(in)
%
% like diff() but also across NaNs
%
% MvdM 2015-05-01 initial version

out = nan(size(in));

nonnan_idx = find(~isnan(in));

in_nonnan = in(nonnan_idx);
in_nonnan_diff = diff(in_nonnan);

out_idx = nonnan_idx(2:end);
out(out_idx) = in_nonnan_diff;

out = out(2:end);
