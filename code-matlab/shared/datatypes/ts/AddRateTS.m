function [TS,gap_idx] = AddRateTS(cfg,TS,tvec)
%ADDRATETS Add usr field to TS that contains the number of occurrances per
% unit time (for S, this is firing rate)
%   [TS,gap_idx] = ADDRATETS(cfg,TS,tvec)
%
%   INPUTS
%      cfg:    Struct with field controlling function behavior (see config
%              options)
%      TS:     TS struct containing timestamps such as spike times
%      tvec:   Time vector (can be obtained from LFP or position data)
%
%   OUTPUTS
%      TS:     TS struct with usr rate field added.
%     gap_idx: Indices where gaps were found (is empty if ~cfg.DetectGaps)
%
%   CONFIG OPTIONS
%
%      cfg.FieldName = 'rate';
%          What to call the usr field that will contain rate data.
%
%      cfg.DetectGaps = 0; 
%          If 1, account for large gaps found in the tvec that might make
%          firing rate incorrectly appear lower, if 0 don't. If you don't
%          have massive gaps, or lots of them, this will not impact the
%          rate very much.
%
%      cfg.method = 'mean';
%          How to determine whether there is a gap. Applies if cfg.DetectGaps
%             'raw'    - if the gap is greater than cfg.threshold (in tvec
%                        units)
%             'mean'   - if the gap is greater than cfg.threshold times
%                        above the mean tvec difference
%             'zscore' - if the gap is greater than cfg.threshold standard
%                        deviations above the mean tvec difference
%
%      cfg.threshold = 2;
%          Threshold for identifying a gap. Applies only if cfg.DetectGaps
%
%      cfg.verbose = 1;
%
% aacarey Jan 2016

% Set config defaults
cfg_def.FieldName = 'rate';
cfg_def.DetectGaps = 0;
cfg_def.method = 'mean';
cfg_def.threshold = 2;
cfg_def.verbose = 1;
mfun = mfilename;

cfg = ProcessConfig(cfg_def,cfg,mfun);

% check inputs
if ~CheckTS(TS)
    error('TS data must have been made with the TS constructor. See ts()')
end

if isempty(TS.t)
    error('TS data is empty. Cannot add rate.')
end

if cfg.verbose; disp([mfun,': Adding rate to usr data']); end

% Get total time
ElapsedTime = tvec(end)-tvec(1);
gap_idx = [];

% Handle gaps in time, if requested
if cfg.DetectGaps
    TimeDifference = diff(tvec);
    
    switch cfg.method
        case 'raw'
            gap_idx = find(TimeDifference >= cfg.threshold);
            
        case 'mean'
            gap_idx = find(TimeDifference >= mean(TimeDifference)*cfg.threshold);
            
        case 'zscore'
            gap_idx = find(zscore(TimeDifference) >= cfg.threshold);
            
        otherwise
            error('Unrecognized option specified in cfg.method')
    end
    
    LostTime = sum(tvec(gap_idx + 1) - tvec(gap_idx));
    ElapsedTime = ElapsedTime - LostTime;     
end

% Calculate rate and add usr field
for iTS = 1:length(TS.t)
    nOccurrances = length(TS.t{iTS});
    TS.usr.(cfg.FieldName)(iTS) = nOccurrances/ElapsedTime;
end

% Record history
TS = History(TS,cfg,mfun);

end

