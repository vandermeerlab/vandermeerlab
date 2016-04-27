function iv_out = InvertIV(iv_in,t0,t1)
% function iv_out = InvertIV(iv_in,t0,t1)
% 
% "flips" included intervals to non-included, and vice versa
%
% MvdM 2016-04-27

iv_out = iv;
iv_out.tstart = cat(1,t0,iv_in.tend);
iv_out.tend = cat(1,iv_in.tstart,t1);

