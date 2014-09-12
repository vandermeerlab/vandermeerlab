% index =   nearest_idx3(values, lookupvalues)
% index =   nearest_idx3(values, lookupvalues, direction)
%
% This function will return the indices of lookupvalues that are closest to
% the numbers in values. 
% NOTE: The function assumes that LOOKUPVALUES are in
% order from smallest to largest. NO WARNING IS given if this is not the
% case!
%
% direction of -1 forces lookup to only look backwards in lookupvalues. If
% direction is 0 (default), however, it chooses the closest indeces. If it
% is 1 it looks foreward in lookupvalues.  If a value falls outside the
% range of lookupvalues, it always chooses the closest index, regardless of
% the value of direction.
%
% for unsorted values, use nearest_idx2.m
%
% MvdM 2014-07-04, based on lookup.c by Mattias Karlsson