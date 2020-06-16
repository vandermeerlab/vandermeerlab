function iv_out = ResizeIVfromTSD(cfg_in, iv_in, tsd_in)
% function iv_out = ResizeIVfromTSD(cfg_in, iv_in, tsd_in)
%
% ResizeIVfromTSD expand or contract intervals to specified threshold of input tsd
% 
% Available modes (set in cfg_in.mode):
%
% * 'lower_threshold': expand each iv to the time at which a LOWER
% cfg_in.threshold than that used for the original iv is crossed.
%
% Example:
%
% original iv (2):    xxxx           xxxx
%             tsd: 001223211000121000222211210
%      new iv (1):   xxxxxxx         xxxxxxxx
%
% NOTE that this mode currently only works if the original ivs were 
% detected with the '>' or '>=' operations.

cfg_def = [];
cfg_def.mode = 'lower_threshold';
cfg_def.threshold = [];

cfg = ProcessConfig(cfg_def, cfg_in);

if isempty(cfg.threshold)
    error('cfg.threshold input is REQUIRED.');
end

iv_out = iv_in;

detec = tsd_in.data > cfg.threshold;
detec = cat(2,0,detec,0); % pad the detection so we can deal with cases where first or last samples are detects

dfs = diff(detec);
up_idx = find(dfs == 1);
down_idx = find(dfs == -1) - 1;

up_t = tsd_in.tvec(up_idx);
down_t = tsd_in.tvec(down_idx);


% loop over ivs to find nearest up- and downcrossings
for iI = length(iv_out.tstart):-1:1
    
    this_center = (iv_out.tstart(iI) + iv_out.tend(iI)) / 2;
    
    up_idx = nearest_idx3(this_center, up_t, -1);
    iv_out.tstart(iI) = up_t(up_idx);
    
    down_idx = nearest_idx3(this_center, down_t, 1);
    iv_out.tend(iI) = down_t(down_idx);
      
end % of loop over ivs


% housekeeping
iv_out.cfg.history.cfg = cat(1, iv_out.cfg.history.cfg,{cfg});
