function c = censored( tau_c, M, T )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% censored - returns the expected false negative rate due to censoring
%
% Usage:
%           c = censored( tau_c, M, T )
%
% Description:  
%   Estimates the percent of a cluster that is lost due to censoring.  Every
% detected spike event turns off spike detection for a period tau_C.  This
% simple function estimates the percent of false negative errors in a cluster
% by calculating what fraction of the entire data set was censored by
% detected events that are external to this cluster.
%
% If tau_c is the duration of the censored period, M is the number of other
% detected events, and T is the duration of the recording period, then the
% fraction of the experiment that was censored is simply tau_c * M / t.
%
% Input: 
%       tau_c -  censored period (ms)
%           M -  number of events detected outside of cluster of interest
%           T -  total duration of experiment (s)
%
% Output:
%       c    - estimated false negative fraction from censoring
%

    c = (tau_c/1000) * M / T;

