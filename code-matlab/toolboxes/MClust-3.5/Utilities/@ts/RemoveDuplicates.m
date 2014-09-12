function [TS_out,nRemoved] = RemoveDuplicates(TS_in)

% @ts/RemoveDuplicates:
% [TS_out,nRemoved] = RemoveDuplicates(TS_in)
%
% MvdM 07

TS_out = TS_in;

n_in = length(TS_in.t);

TS_out.t = unique(TS_in.t);
n_out = length(TS_out.t);

nRemoved = n_in - n_out;



