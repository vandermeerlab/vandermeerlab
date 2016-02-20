function [pass_flag,badFields] = CheckConfig(cfg_def,cfg_in,verbose,varargin)
%CHECKCONFIG Identify unrecognized cfg_in fields
%   pass_flag = CheckConfig(cfg_def,cfg_in,verbose) displays fields in 
%   cfg_in that are not fields in cfg_def. Returns a pass flag that can be 
%   used to hault further analyses.
% 
%   pass_flag = CheckConfig(cfg_def,cfg_in,verbose,callername) if verbose =
%   1 and callername is present (name of invoking function), then the
%   warning message will be more informative.
%  
%   Why use this function?
%   If you have a ocnfig feild that is speled incorreclty, the function
%   or script will use the default config field [potentially] without you ever 
%   finding out! Using a check will help avoid this situation. 
%   Alternatively, it may point out that you have used a config intended
%   for another function for the current function (maybe you forget to
%   reset your config).
%
%   INPUTS:
%   cfg_def: default config containing all the fields the function may use
%   cfg_in:  config to be checked 
%   verbose: if 1, display incorrect fields; if 0, don't
%   callername: the name of the invoking function (see mfilename)
%
%   OUTPUTS:
%   pass_flag: 1 if all cfg_in fields are also cfg_def fields, 0 otherwise
%   badFields: cell array of unrecognized fields
%
% see also ProcessConfig 
%
%   aacarey Oct 2015

%%

if isempty(cfg_in)
    cfg_in = cfg_def;
end

pass_flag = 1;

if ~isempty(cfg_in)
    inFields = fieldnames(cfg_in);
    badFields = inFields(arrayfun(@(x) ~isfield(cfg_def,x),inFields));
    
    if ~isempty(badFields) && verbose
        in_mfun = '';
        if ~isempty(varargin) && ischar(varargin{1})
            in_mfun = [' in ',varargin{1}];
        end
        disp(['WARNING',in_mfun,': unrecognized config field(s) ',strjoin(badFields',', ')])
    end
end

end

