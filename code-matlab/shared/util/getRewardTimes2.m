function reward_t = getRewardTimes2()
% function reward_t = getRewardTimes()
%
% load events file and extract times of first feeder fire
%
% tested on Multiple-T and Multiple-T ARL data (Redish lab)
%
% must be run within session folder

LoadExpKeys;
evt = LoadEvents([]);

keep = ~cellfun('isempty',evt.label); evt = SelectTS([],evt,keep);

if isfield(ExpKeys,'FeederL1') % Multiple-T ARL
    
    feeders = cat(2, ExpKeys.FeederL1, ExpKeys.FeederR1);
    reward_t = [];
    ll = @(x) x(end); % function to get last character of input
    for iF = 1:length(feeders)
        
        keep_idx = find(num2str(feeders(iF)) == cellfun(ll, evt.label));
        reward_t = cat(1, reward_t, evt.t{keep_idx});
        
    end
    reward_t = sort(reward_t);
    keep = find(diff(reward_t) > 1); keep = [1; keep + 1];
    reward_t = reward_t(keep);
    
elseif isfield(ExpKeys,'Feeder2') % Multiple-T
    
    reward_t = [];
    ll = @(x) x(end); % function to get last character of input
    keep_idx = find(num2str(ExpKeys.Feeder2) == cellfun(ll, evt.label));
    reward_t = cat(1, reward_t, evt.t{keep_idx});
    
    reward_t = sort(reward_t);
    keep = find(diff(reward_t) > 1); keep = [1; keep + 1];
    reward_t = reward_t(keep);
else
    warning('Could not load reward times.');
    reward_t = [];
end