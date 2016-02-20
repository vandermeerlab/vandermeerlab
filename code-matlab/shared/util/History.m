function var_out = History(var_in,mfun,cfg)
%HISTORY Add or update a variable's config history
%
% var_out = HISTORY(var_in,mfun,cfg)
%
%   INPUTS
%      var_in: the variable you want to add or update history for
%      mfun: name of caller function
%      cfg: the config you want to store in variable history
%
%   OUTPUT
%      var_out: variable with updated config history in cfg.history field:
%             .cfg.history.mfun = nx1 cell array of strings containing
%                     previous function names. {1,1} is the first function,
%                     {2,1} is the second function, and so on
%             .cfg.history.cfg = nx1 cell array of structs containing
%                     previous configs corresponding to the entries in the
%                     mfun field. {1,1} is the config for the first function,
%                     {2,1} is for the second function, and so on
%
% see also GetHistory
%
% Dec 2015

var_out = var_in;

if isfield(var_out,'cfg') && isfield(var_out.cfg,'history')
    var_out.cfg.history.mfun = cat(1,var_out.cfg.history.mfun,mfun);
    var_out.cfg.history.cfg = cat(1,var_out.cfg.history.cfg,{cfg});
else
    var_out.cfg.history.mfun{1} = mfilename;
    var_out.cfg.history.cfg{1} = [];
end
end

