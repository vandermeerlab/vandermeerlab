function [iva,idxa] = OverlapIV(cfg_in,ivA,ivB)
%OVERLAPIV Return the intervals and indices of ivA that span the same 
% regions as the intervals in ivB.
%
%   [iva,idxa] = OVERLAPIV(cfg,ivA,ivB)
%                                                         
%   ivA        ___________________  _  __________       _____   _
%
%   ivB     _________       ____    _____          __         ______ ___
%
%   iva        ___________________  _  __________               _
%
%   idxa         1            1     2   3                       5 
%
%    INPUTS
%       cfg: config struct with fields controlling function behaviour
%       ivA: interval data to be compared to ivB
%       ivB: interval data for comparison
%   
%    OUTPUTS
%       iva: interval data from ivA that overlaps with the intervals in ivB
%      idxa: the indices (identities) of the overlapping intervals. Note
%            that the number of indices returned will be the same regardless 
%            of the order of the input ivs. So if you reverse the input order
%            of ivA and ivB, then the indices returned will tell you the 
%            identities of the intervals that have been paired up with each 
%            other. In the ASCII example above, the indices returned are
%            [1,1,2,3,5] if ivA is the first input; if instead ivB is the
%            first input, the indices would be [1,2,3,3,5]. In this way you
%            know that the first interval in A overlaps with both the first
%            and second intervals in B.
%
%    CONFIG OPTIONS
%       cfg.verbose = 1; If 1, tell me how many intervals came in and how
%            many went out; if 0, don't
%
% aacarey Nov 2015

%%

mfun = mfilename;

% initialize default parameters
cfg_def.verbose = 1;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

idxa = [nan(size(ivA.tstart)); nan(size(ivB.tstart))]; % allocate space
current_idx = 1;
for iB = 1:length(ivB.tstart)
    for iA = 1:length(ivA.tstart)
        if (ivB.tstart(iB) >= ivA.tstart(iA) && ivB.tstart(iB) <= ivA.tend(iA))...
                || (ivB.tend(iB) >= ivA.tstart(iA) && ivB.tend(iB) <= ivA.tend(iA))...
                || (ivA.tstart(iA) >= ivB.tstart(iB) && ivA.tstart(iA) <= ivB.tend(iB))...
                || (ivA.tend(iA) >= ivB.tstart(iB) && ivA.tend(iA) <= ivB.tend(iB))
            
            idxa(current_idx) = iA;
            current_idx = current_idx +1;
        end
    end
end
idxa = idxa(~isnan(idxa));

% select intervals to keep
cfg_temp = []; cfg_temp.verbose = 0;
iva = SelectIV(cfg_temp,ivA,idxa);

% remove doubles (this happens when more than one ivB interval corresponds to the same ivA interval)
%   ivA   ____________
%   ivB    ____   ____

cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.rmdoubles = 1;
iva = RemoveIV(cfg_temp,iva);

% talk to me
if cfg.verbose
    disp([mfun,': ',num2str(length(ivA.tstart)),' intervals in, ',num2str(length(iva.tstart)),' intervals out.'])
end

% housekeeping
iva.cfg.history = ivA.cfg.history;
iva = History(iva,mfun,cfg);

end
