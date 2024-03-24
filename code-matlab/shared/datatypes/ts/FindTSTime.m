function [t,idx] = FindTSTime(cfg_in,ts_in,t_in)
% function [t,idx] = FindTSTime(cfg,ts_in,t_in)
%
% returns the timestamp and corresponding index of 
% either the next or previous time, across all ts_in.t{} 
% timestamps given some time t_in
%
% cfg_in.mode = 'next'; % 'prev'
% 
% see also FindFieldTime()
%
% MvdM 2023, based on FindFieldTime()

if length(ts_in) ~= 1
   error('ts_in must be length 1 (is %d)',length(ts_in)); 
end

cfg = [];
cfg.mode = 'next'; % 'prev'
cfg = ProcessConfig(cfg,cfg_in);

% assemble arrays
all_label_idxs = [];
all_times = [];

for iF = 1:length(ts_in.t)
    
    curr_t = ts_in.t{iF};
    curr_t = curr_t(:); % convert to column vector
    
    all_times = cat(1,all_times,curr_t);
    all_label_idxs = cat(1,all_label_idxs,iF.*ones(size(curr_t)));
    
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
all_label_idxs = all_label_idxs(keep_idx);

% check if still times remaining
if isempty(all_times)
   t = []; idx = [];
   warning('No events found.');
   return;
end

% sort
[all_times,sort_idx] = sort(all_times,'ascend');
all_label_idxs = all_label_idxs(sort_idx);

% retrieve target item
switch cfg.mode
    case 'next'
        t = all_times(1);
        idx = all_label_idxs(1);
    case 'prev'
        t = all_times(end);
        idx = all_label_idxs(end);     
end
