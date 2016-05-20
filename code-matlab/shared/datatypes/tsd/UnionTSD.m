function tsd_out = UnionTSD(cfg_in,tsd1,tsd2)
% function tsd_out = UnionTSD(cfg_in,tsd1,tsd2)
%
% returns the union (merge) of two TSD inputs
%
% example:
%
% tsd1.tvec = [1 3 5]; tsd1.data = [1 3 5]; 
%
% tsd2.tvec = [2 4]; tsd1.data = [22 44];
%
% tsd_out.tvec = [1 2 3 4 5]; tsd_out.data = [1 22 3 44 5];
%
% - throws error if any of the tvec values are identical
% - output tvec will be sorted
%
% MvdM 2016-04-15 initial version
% NOTE merging of usr fields doesn't work well if one input is empty

cfg_def = [];
cfg = ProcessConfig2(cfg_def,cfg_in);

mfun = mfilename;

%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

if ~CheckTSD(tsd1) | ~CheckTSD(tsd2)
   error('Malformed TSD'); 
end

% check if any timestamps are identical; if so, give up
common_values = intersect(tsd1.tvec,tsd2.tvec);
if ~isempty(common_values)
   error('tsd1 and tsd2 have common tvec value(s)'); 
end

% check if any usr field names are different; if so, give up
fn1 = []; fn2 = [];
if xor(isfield(tsd1,'usr'),isfield(tsd2,'usr')) & (~isempty(tsd1.tvec) & ~isempty(tsd2.tvec)) % OK if one is empty...
    error('tsd1/tsd2 usr field presence mismatch');
elseif isfield(tsd1,'usr') & isfield(tsd2,'usr')
    fn1 = fieldnames(tsd1.usr); fn2 = fieldnames(tsd2.usr);
    fnx = setxor(fn1,fn2); % names that are in one, but not the other
    if ~isempty(fnx)
        error('tsd1 and tsd2 have different usr field names');
    end
end

%%%%%%%%%%%%%
%%% MERGE %%%
%%%%%%%%%%%%%

tsd_out = tsd; % create blank tsd

% construct new row tvec
if ~isrow(tsd1.tvec), tsd1.tvec = tsd1.tvec'; end
if ~isrow(tsd2.tvec), tsd2.tvec = tsd2.tvec'; end

tsd_out.tvec = cat(2,tsd1.tvec,tsd2.tvec);

% construct new data -- first check shape
if size(tsd1.data,1) > size(tsd1.data,2)
    fprintf('WARNING: UnionTSD.m: tsd1 has more nSignals than nSamples!\n');
end

if size(tsd2.data,1) > size(tsd2.data,2)
    fprintf('WARNING: UnionTSD.m: tsd2 has more nSignals than nSamples!\n')
end

tsd_out.data = cat(2,tsd1.data,tsd2.data);

[tsd_out.tvec,sort_idx] = sort(tsd_out.tvec,'ascend');
tsd_out.data = tsd_out.data(:,sort_idx);

% handle usr fields
for iF = 1:length(fn1)
   
    tsd_out.usr.(fn1{iF}) = cat(2,tsd1.usr.(fn1{iF}),tsd2.usr.(fn1{iF}));
    tsd_out.usr.(fn1{iF}) = tsd_out.usr.(fn1{iF})(:,sort_idx);
    
end

% housekeeping
tsd_out.cfg.history.mfun = cat(1,tsd_out.cfg.history.mfun,mfun);
tsd_out.cfg.history.cfg = cat(1,tsd_out.cfg.history.cfg,{cfg});

