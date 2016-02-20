function epochTimes = getEpochTimes(cfg_in,ExpKeys,metadata)
% function epochTimes = getEpochTimes(cfg,ExpKeys,metadata)
%
% returns iv cell array defining task epochs for Alyssa's T-maze task, such
% as first 5 trials, etc...
%
% INPUTS:
% cfg_def.labels = {}; % e.g. {'pre','pre_start','pre_end','taskrest_start','taskrest_end','post'};
% cfg_def.nTrials = 5; % number of trials to use for taskrest* epochs
%
% OUTPUTS:
%
% epochTimes: 1 x nLabels cell array, each cell contains iv with start and
%  end times for that epoch
%
% EXAMPLE USAGE:
% cfg = []; cfg.labels = {'taskrest_start'};
% getEpochTimes(cfg,ExpKeys,metadata)
%
% MvdM 2015-10-01 initial version

cfg_def = {};
cfg_def.nTrials = 5;

cfg = ProcessConfig2(cfg_def,cfg_in);

epochTimes = {};
for iL = length(cfg.labels):-1:1
    
    switch cfg.labels{iL}
        
        case 'pre' % full prerecord
            
            epochTimes{iL} = iv(ExpKeys.prerecord);
            
        case 'pre_start' % first half of prerecord
            
            epochTimes{iL} = iv([ExpKeys.prerecord(1) ExpKeys.prerecord(1)+(ExpKeys.prerecord(2)-ExpKeys.prerecord(1))/2]);
            
        case 'pre_end' % second half of prerecord
    
            epochTimes{iL} = iv([ExpKeys.prerecord(1)+(ExpKeys.prerecord(2)-ExpKeys.prerecord(1))/2 ExpKeys.prerecord(2)]);
            
        case 'post' % full postrecord
    
            epochTimes{iL} = iv(ExpKeys.postrecord);
            
        case 'taskrest_start' % first cfg.nTrials intertrial rest periods
            
            nTrialsTotal = length(metadata.taskvars.rest_iv.tstart);
            temp_end = min(nTrialsTotal,cfg.nTrials);
            
            epochTimes{iL} = iv(metadata.taskvars.rest_iv.tstart(1:temp_end),metadata.taskvars.rest_iv.tend(1:temp_end));
            
        case 'taskrest_end' % last cfg.nTrials intertrial rest periods
    
            nTrialsTotal = length(metadata.taskvars.rest_iv.tstart);
            temp_start = max(1,nTrialsTotal-cfg.nTrials+1);
            
            epochTimes{iL} = iv(metadata.taskvars.rest_iv.tstart(temp_start:end),metadata.taskvars.rest_iv.tend(temp_start:end));
    end
end
