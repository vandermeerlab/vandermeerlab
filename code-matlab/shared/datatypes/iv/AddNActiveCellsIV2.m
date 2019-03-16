function iv_in = AddNActiveCellsIV2(cfg_in,iv_in,Q)
% function iv_out = AddNActiveCellsIV(cfg_in,iv_in,S,Q)
%
% Use this version if you need to repeatedly call AddNACtive CellsIV in a
% loop. MakeQfromS is a slow function when computed over many shuffles.
% If S is identical across all shuffles, you can compute Q a single time
% outside of the loop and pass it into MakeQfromS2. This substantially
% decreases the redundancy.
%
%   INPUT:
%       cfg_in - input cfg parameters
%       iv_in - input iv
%       Q - Q matrix (see MakeQfromS)
%
%   OUTPUT:
%       iv_out - iv with nActiveCells usr field added
%
%   CFG OPTIONS:
%
%       cfg_def.label = 'nActiveCells'; What to call the field that
%                    contains this data.
%       cfg_def.dt = 0.001; % bin size for binning spikes
%
% mvdm 2015-02-01 initial version
% note: uses binning method for speed, i.e. is not exact

% parse cfg parameters
cfg_def = [];
cfg_def.label = 'nActiveCells';
cfg_def.dt = 0.005;
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);


% check inputs
if ~CheckIV(iv_in)
    error('Interval data must have been made with the iv constructor')
end

if isempty(iv_in.tstart)
    error('Interval data is empty')
end

if cfg.verbose; disp([mfun,': Adding the number of active units during each interval']); end

% make matrix of spike counts (nCells x nTimeBins) 
%cfg_temp.dt = cfg.dt;
%cfg_temp.verbose = 0;
%Q = MakeQfromS(cfg_temp,S);

% get number of active cells for each
for iIV = length(iv_in.tstart):-1:1
    
    Qr = restrict(Q,iv_in.tstart(iIV),iv_in.tend(iIV));
    
    if ~isempty(Qr.tvec)
        spk_counts = sum(Qr.data,2);
        nC(iIV) = sum(spk_counts > 0);
    else % no spikes in this iv
        nC(iIV) = 0;
    end
    
end

if ~isfield(iv_in,'usr') 
    iv_in.usr = [];
end

if cfg.verbose && isfield(iv_in.usr,cfg.label)
    disp(['WARNING in ',mfun,': iv data already includes usr.',cfg.label,', overwriting...'])
end
iv_in.usr.(cfg.label) = nC(:); % ensure column shape


% housekeeping
iv_in.cfg.history.mfun = cat(1,iv_in.cfg.history.mfun,mfun);
iv_in.cfg.history.cfg = cat(1,iv_in.cfg.history.cfg,{cfg});