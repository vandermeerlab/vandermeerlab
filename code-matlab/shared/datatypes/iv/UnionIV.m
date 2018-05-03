function iv1 = UnionIV(cfg_in,iv1,iv2)
% function iv = UnionIV(cfg,iv1,iv2)
%
% union of iv objects
%
% output is resorted based on start times (ascending)
%
% MvdM 2014-08-28 initial version

cfg_def = [];
cfg = ProcessConfig(cfg_def,cfg_in);

mfun = mfilename;

if isempty(iv1) | isempty(iv1.tstart) % function should work for empty arguments
    %iv1 = iv([],[]);
    iv1 = iv2; return;
end

if isempty(iv2) | isempty(iv2.tstart)
    %iv2 = iv([],[]);
    return;
end

if ~CheckIV(iv1) | ~CheckIV(iv2)
   error('Malformed IV.'); 
end

iv1.tstart = cat(1,iv1.tstart,iv2.tstart); % note that constructor guarantees column vectors
iv1.tend = cat(1,iv1.tend,iv2.tend);

[iv1.tstart,sort_idx] = sort(iv1.tstart,'ascend');
iv1.tend = iv1.tend(sort_idx);

% process usr fields -- needs check for columnness, N-D cases etc..
if (isfield(iv1,'usr') && ~isempty(iv1.usr)) && (isfield(iv2,'usr') && ~isempty(iv2.usr))
    ivfields = fieldnames(iv1.usr);
    for iField = 1:length(ivfields)
        iv1.usr.(ivfields{iField}) = cat(1,iv1.usr.(ivfields{iField}),iv2.usr.(ivfields{iField})); % note, assumes column shape!
        iv1.usr.(ivfields{iField}) = iv1.usr.(ivfields{iField})(sort_idx);
    end
elseif (~isfield(iv1,'usr') && isfield(iv2,'usr')) | (isfield(iv1,'usr') && ~isfield(iv2,'usr'))
    error('Cannot merge usr fields: only one input has an usr field');
elseif (isfield(iv1,'usr') && isempty(iv1.usr)) && (isfield(iv2,'usr') && ~isempty(iv2.usr))
    error('Cannot merge usr fields: only iv2 has non-empty usr field');
elseif (isfield(iv1,'usr') && ~isempty(iv1.usr)) && (isfield(iv2,'usr') && isempty(iv2.usr))
    error('Cannot merge usr fields: only iv1 has non-empty usr field');
end

% housekeeping
iv1.cfg.history.mfun = cat(1,iv1.cfg.history.mfun,mfun);
iv1.cfg.history.cfg = cat(1,iv1.cfg.history.cfg,{cfg});