function ts_out = concatenateTS(ts1,ts2)
%concatenateTS concatenate two ts objects
%   Concatenates .t fields and all related information (labels, usr fields)
%
%   ts.out = concatenateTS(ts1,ts2)
%
%   INPUTS
%      ts1  - first ts (second ts will be appended to this)
%      ts2  - second ts
%
%   OUTPUTS
%      ts_out - concatenated ts
%
% MvdM 2015-11-18 initial version

mfun = mfilename; % get name of function as string

% append t and label
ts_out = ts1;

ts_out.t = cat(2,ts_out.t,ts2.t);
ts_out.label = cat(2,ts_out.label,ts2.label);

% prepare to append data from usr fields

% first check if have same usr fields
if ~isfield(ts_out,'usr')
    nFields1 = 0;
else
    ts1_fieldnames = fieldnames(ts_out.usr);
    nFields1 = length(ts1_fieldnames);
end

if ~isfield(ts2,'usr')
    nFields2 = 0;
else
    ts2_fieldnames = fieldnames(ts2.usr);
    nFields2 = length(ts2_fieldnames);
end
    
if nFields1 ~= nFields2
    error('Inputs have non-matching usr field counts.');
elseif sum(strcmp(sort(ts1_fieldnames),sort(ts2_fieldnames))) ~= nFields1 % some name mismatch
    error('Inputs have non-matching usr field names.');
end

% append
if isfield(ts_out,'usr') && ~isempty(ts_out.usr)
    ts1_fieldnames = fieldnames(ts_out.usr);
    for iField = 1:length(ts1_fieldnames)
        ts_out.usr.(ts1_fieldnames{iField}) = cat(2,ts_out.usr.(ts1_fieldnames{iField}),ts2.usr.(ts1_fieldnames{iField}));
    end
end

% keep a record of cfg history -- note this is a case with two parents, not
% handled well
if isfield(ts_out,'cfg')
    ts_out.cfg.history.mfun = cat(1,ts_out.cfg.history.mfun,mfilename);
    ts_out.cfg.history.cfg = cat(1,ts_out.cfg.history.cfg,{});
else
    ts_out.cfg.history.mfun = mfun;
    ts_out.cfg.history.cfg = {};
end


