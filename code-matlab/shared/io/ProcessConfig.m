function [cfg_out,pass_flag,badFields] = ProcessConfig(cfg_def,cfg_in,varargin)
%% PROCESSCONFIG Replace default parameters in cfg_def with parameters in 
% cfg_in and spell-check incoming config fields.
% A config, typically named cfg, is a struct containing fields that control
% function behavior. Functions that use configs should initialize all config 
% fields (cfg_def) within the function body. 
% For information on how to initialize a config inside of your function, 
% open ProcessConfig and read the "Help" section.
%
%  [cfg_out] = PROCESSCONFIG(cfg_def,cfg_in,)
%  [cfg_out] = PROCESSCONFIG(cfg_def,cfg_in,callername)
%
%   INPUTS:
%    required:
%       cfg_def: default parameters stored as a struct within the caller function
%       cfg_in: input parameters from the input cfg struct
%    varargins:
%       callername: the name of the invoking function (see mfilename) for
%         more informative spell-checking
%
%   OUTPUT:
%       cfg_out : Returns a new cfg struct containing the default parameters and any updated cfg
%       parameters given the input struct.
%
%    additional (optional) output parameters:
%     [cfg_out,pass_flag,badFields] = ProcessConfig(cfg_def,cfg_in,varargin)
%      pass_flag: 1 if all cfg_in fields are also cfg_def fields, 0 otherwise
%      badFields: cell array of unrecognized fields
%
%   Recommended: initialize cfg_def.verbose = 1 inside your function. Why?
%   If you have a ocnfig feild that is speled incorreclty, the function
%   or script will use the default config field [potentially] without you ever 
%   finding out! If you allow ProcessConfig to be verbose, it will notify you if
%   it does not recognize a field (for example, if the field is actually called
%   'threshold' and you typed 'thershold'). 
%   In some cases, it may help point out that you have used a config intended
%   for another function (maybe you forget to reset your config and it just 
%   happens to have an unrecognized field).
%
%   see also ProcessConfig2, CheckConfig
% 
% youkitan 2014-11-04 (inital version ProcessConfig2)
% aacarey Nov 2015 (config spell checker, more info, renamed to ProcessConfig)

%% Help

% ~~~ How to set config defaults within your function ~~~~~~~~~~~~~~~~~~~~~
%
% function out = myFun(cfg_in,Arg1,Arg2)
%
% mfun = mfilename; % this grabs the name of the caller function
%
% % Initialize default parameters here (make sure you give them more useful
% names than 'parameter1')!
% cfg_def.verbose = 1; % by default, display helpful text 
% cfg_def.parameter1 = 5; 
% cfg_def.parameter2 = []; 
% cfg_def.parameter3 = 'interpolate';
%
% % Now reassign config parameters according to cfg_in
% cfg = ProcessConfig(cfg_def,cfg_in,mfun);
% 
% % do function stuff here
% if cfg.verbose
%   disp('It was a sunny but cold day at Dartmouth when this was written')
% end
% % more stuff here
% 
% end

% ~~~ The importance of refreshing your config between functions ~~~~~~~~~~ 
%
% Sometimes functions use the same fieldnames. Consider LoadSpikes() and
% LoadCSC(): They both have fieldnames cfg.fc specifying filenames to load.
% Early on when the codebase was still young, LoadSpikes() did not perform 
% a check on the  filename it tried to load. So, if you did not set cfg = [] 
% and then reset cfg.fc between LoadCSC and LoadSpikes, LoadSpikes
% would try to load a .ncs file (in fact, it would try to load just about
% anything you put in). This has been fixed, but there are still other functions 
% that share config fieldnames. For example, many functions have cfg.threshold. If you
% wanted to use cfg_def.threshold = 5 in function2(), but you had set
% cfg.threshold = 1 in function1(), and then you  passed function1's config in
% function2, the results would not be as you expect. (And that previous
% sentence is why scientists don't teach English classes).
% TLDR: always set cfg = [] between function calls (or use different config
% names).

%% Set output cfg to default cfg

cfg_out = cfg_def;

%% Process input cfg parameters into workspace
% NOTE - The two processes (replace and add) are separated only for clarity 
% and ease of debugging. It can easily be compressed to half the amount of code (and
% therefore time).

%~~~~~~~~~~ PROCESS CONFIG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ~isempty(cfg_in)
    
    % If there are default parameters,replace them with cfg_in input
    if ~isempty(cfg_def)
        default_F = fieldnames(cfg_def);
        numF = length(default_F);

        for i = 1:numF
            idF = default_F{i};
            if isfield(cfg_in,idF)
                cfg_out.(idF) = cfg_in.(idF);
            end %set cfg_out fields    
        end %iterate default cfg fields
    end

    % Add input parameters
    new_F = fieldnames(cfg_in);
    numF = length(new_F);
    
    for i = 1:numF
        inF = new_F{i};
        if ~isfield(cfg_def,inF)
            cfg_out.(inF) = cfg_in.(inF);
        end %set cfg_out fields    
    end %iterate extra cfg fields
 
%~~~~~~~~~~ CHECK CONFIG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (isfield(cfg_out,'verbose') && cfg_out.verbose) || nargout > 2
    if isempty(cfg_in)
        cfg_in = cfg_def;
    end
    
    pass_flag = 1;
    
    if ~isempty(cfg_in)
        inFields = fieldnames(cfg_in);
        badFields = inFields(arrayfun(@(x) ~isfield(cfg_def,x),inFields));
        
        if ~isempty(badFields) && isfield(cfg_out,'verbose') && cfg_out.verbose
            
            in_mfun = ''; % set empty in case callername isn't given
            
            if ~isempty(varargin) && ischar(varargin{1})
                in_mfun = [' in ',varargin{1}];
            end
            disp(['WARNING',in_mfun,': unrecognized config field(s) ',strjoin(badFields',', ')])
        end
    end
end
end
