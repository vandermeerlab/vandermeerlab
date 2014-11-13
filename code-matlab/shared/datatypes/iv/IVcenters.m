function c = IVcenters(iv_in)
% returns nIV x 1 vector with center times of iv_in

if ~CheckIV(iv_in)
    return;
end

c = nanmean(cat(2,iv_in.tstart,iv_in.tend),2);