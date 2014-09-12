% ProcessConfig.m
% 
% NOT A FUNCTION -- allows for access to current workspace
%
% MvdM 2014-06-17
% modified 2014-06-25 to separate cfg_in from cfg

if isempty(cfg_in)
    return;
end

TEMP_f = fieldnames(cfg_in);
nF = length(TEMP_f);

for iF  = 1:nF
   
    cur_f = TEMP_f{iF};
    eval(sprintf('cfg.%s = cfg_in.%s;',cur_f,cur_f));
    
end
clear TEMP_F