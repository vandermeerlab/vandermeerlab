function [history,idx] = GetHistory(cfg_in,var_in)
%GETHISTORY Extract config history from datatype variable
%
% history = GETHISTORY(cfg,var_in)
%
%   INPUTS
%     cfg:      config struct with fields controlling THIS function's behavior
%     var_in:   datatype variable with config history you want to extract
%
%    OUTPUTS
%      history: the config corresponding to the target function name, or
%               a specific parameter requested in cfg.parameter (see config
%               options).
%      idx:     the indices of the targets in config history. If cfg.target
%               = 'all', then idx = [].
%
%    CONFIG OPTIONS
%     cfg.target = 'all'; function name specifying which config to grab
%               'all'  -  return entire history. In this case history has
%                         two fields: history.mfun and history.cfg
%               mfun   -  name of specific function (ex: 'TSDtoIV',
%                         'FilterLFP') for the config you want returned. In
%                         this case, history is the exact config requested.
%     cfg.parameter = ''; String specifying single parameter in target
%                         config to return (ex: if yout put in 'threshold' 
%                         returns threshold from target config instead of 
%                         entire config)
%     cfg.verbose = 1;
%
% see also History, rmHistory
%
% aacarey Oct 2015, edit Dec 2015

cfg_def.target = 'all';
cfg_def.parameter = '';
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

switch cfg.target
    case 'all'
        history = var_in.cfg.history;
        if nargout == 2
            idx = [];
        end
    otherwise
        % get the names of all the functions in the history
        mfun_names = var_in.cfg.history.mfun;
        
        % get all the configs in the history
        mfun_cfgs = var_in.cfg.history.cfg;
        
        % select the target config
        keep = strcmp(cfg.target,mfun_names);
        if cfg.verbose
            fprintf('%s: Target function ''%s'' has %d call(s) recorded in config history\n',mfun,cfg.target,sum(keep));
        end
        history = mfun_cfgs(keep);
        
        if nargout == 2
            idx = find(keep);
        end
        
        % if requested, return only the value(s) stored in a specific target config
        % field
        if isempty(history) && ~isempty(cfg.parameter)
            error('Target config for %s is empty: cannot return parameter ''%s''',cfg.target,cfg.parameter)
        elseif ~isempty(history) && ~isempty(cfg.parameter)
            for iCFG = 1:length(history)
                history{iCFG} = history{iCFG}.(cfg.parameter);
            end
        end
end

end

