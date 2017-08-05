function [t,fieldname] = FindFieldTime(cfg_in,t_struct,t_in)
% function [t,fieldname] = FindFieldTime(cfg,t_struct,t_in)
%
% assuming a t_struct with fields containing event timestamps (e.g.
% x.reward_time, x.trial_start) returns the timestamp and corresponding field
% name of either the next or previous event given some timestamp t_in
%
% cfg_in.mode = 'next'; % 'prev'
%
% MvdM 2015

if length(t_in) ~= 1
   error('t_in must be length 1 (is %d)',length(t_in)); 
end

cfg = [];
cfg.mode = 'next'; % 'prev'
cfg = ProcessConfig(cfg,cfg_in);

if ~isfield(cfg,'fields')
    cfg.fields = fieldnames(t_struct);
end

% assemble arrays
all_fieldIDs = [];
all_times = [];

for iF = 1:length(cfg.fields)
    
    curr_t = getfield(t_struct,cfg.fields{iF});
    curr_t = curr_t(:); % convert to column vector
    
    all_times = cat(1,all_times,curr_t);
    all_fieldIDs = cat(1,all_fieldIDs,iF.*ones(size(curr_t)));
    
end

% keep only required items
switch cfg.mode
    case 'next'
        keep_idx = find(all_times > t_in);
    case 'prev'
        keep_idx = find(all_times < t_in);
    otherwise
        error('Unknown mode %s',cfg.mode)
end
all_times = all_times(keep_idx);
all_fieldIDs = all_fieldIDs(keep_idx);

% check if still times remaining
if isempty(all_times)
   t = []; fieldname = {};
   return;
end

% sort
[all_times,sort_idx] = sort(all_times,'ascend');
all_fieldIDs = all_fieldIDs(sort_idx);

% retrieve target item
switch cfg.mode
    case 'next'
        t = all_times(1);
        fieldname = cfg.fields(all_fieldIDs(1));
    case 'prev'
        t = all_times(end);
        fieldname = cfg.fields(all_fieldIDs(end));     
end
