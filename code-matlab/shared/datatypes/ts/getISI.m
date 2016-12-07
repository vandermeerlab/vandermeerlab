function isi = getISI(cfg_in,S)
%% GETISI Calculate interspike intervals
%   isi = getISI(cfg_in,S) returns the ISI values for each cell in a TS variable
%
%   INPUTS:
%       cfg_in: config field
%            S: ts struct
%
%   OUTPUTS:
%       isi: {NxI} cell array where N is the number of cells and I is the number of
%            intervals set by the cfg.iv parameter (default 1). Each cell contains the
%            interspike (event) intervals as a [Dx1] double.
%
%   CONFIG OPTIONS
%       cfg.iv = []; If specified, return ISIs for each interval (output is nCells x
%                nIntervals). Otherwise the interval is the whole TS time span.
%
%   See also GETSPIKECOUNT
%
% youkitan 2016-07-28 initial version

%% Parse cfg parameters and error check

cfg_def = [];
cfg_def.iv = []; % if specified, return spike counts for each interval (output is nCells x nIntervals)

cfg = ProcessConfig(cfg_def,cfg_in);

if ~CheckTS(S)
   error('Input is not a correctly formed ts.'); 
end

%% do the things

% create dummy iv to avoid repetition later
if isempty(cfg.iv)
   cfg.iv = iv();
   cfg.iv.tstart = firstSpike(S)-eps;
   cfg.iv.tend = lastSpike(S)+eps;
end

% create array for efficiency
nCells = length(S.t);
nIV = length(cfg.iv.tstart);
isi = cell(nCells,nIV);

% iterate and get ISIs
for iI = 1:nIV % for each interval, get spike counts
   this_S = restrict(S,cfg.iv.tstart(iI),cfg.iv.tend(iI));
   
   for iC = 1:nCells
       isi{iC,iI} = diff(this_S.t{iC});
   end  % of cells
    
end % of intervals
