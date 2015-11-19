function S_out = OrderSelectS(cfg_in,S,idx)
%ORDERSELECTS Order and/or select spiketrains 
%   Return spiketrains and all related information (labels, usr fields)
%   according to the specified indices in idx.
%
%   S = OrderSelectS(cfg,S,idx)
%
%   INPUTS
%      S   - spiketrains (output from LoadSpikes)
%      idx - indices specifying which spiketrains to keep and/or how to
%            redorder them
%            ex: if idx = [54 5 17] then S_out.t(1) corresponds to S.t(54)
%            and so on
%
%   OUTPUTS
%      S_out - spiketrains ordered or selected according to idx
%
%   CONFIG OPTIONS
%      cfg.verbose = 1; % if 1, tell me how many spiketrains went in and
%      how many went out; if 0, don't
%
%   NOTE: this function makes a copy of S and then reorders or selects S.t,
%   S.label, and S.usr.data. It also updates cfg history. This function does 
%   not know about any additional data fields that may be present.
%
% aacarey Oct 2015, edit Nov 2015

%% Parse cfg parameters

cfg_def.verbose = 1;

mfun = mfilename; % get name of function as string
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% do the thing
S_out = S;
S_out.t = S_out.t(idx);
S_out.label = S_out.label(idx);
%S_out.usr.data = S_out.usr.data(idx);

% also select data from other same-length usr fields
if isfield(S_out,'usr') && ~isempty(S_out.usr)
    Sfields = fieldnames(S_out.usr);
    for iField = 1:length(Sfields)
        S_out.usr.(Sfields{iField}) = S_out.usr.(Sfields{iField})(idx);
    end
end

if cfg.verbose % talk to me
    disp([mfun,': ',num2str(length(S.t)),' spiketrains in, ',num2str(length(S_out.t)),' spiketrains out'])
end

% keep a record of cfg history

if isfield(S,'cfg')
    S_out.cfg = S.cfg;
    S_out.cfg.history.mfun = cat(1,S.cfg.history.mfun,mfilename);
    S_out.cfg.history.cfg = cat(1,S.cfg.history.cfg,{cfg});
else
    S_out.cfg.history.mfun = mfun;
    S_out.cfg.history.cfg = {cfg};
end

end

