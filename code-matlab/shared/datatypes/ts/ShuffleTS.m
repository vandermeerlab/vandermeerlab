function ts_in = ShuffleTS(cfg_in,ts_in)
% function ts_out = ShuffleTS(cfg,ts_in)
%
% shuffle spikes in ts struct
%
% INPUTS:
%
% ts_in: input ts
%
% OUTPUT:
%
% ts_out: shuffled ts
%
% CONFIGS:
%
% cfg_def.mode = 1; % 1: shuffle labels, 2: shuffle ISIs, 3: uniform random
%   with same number of spikes
% cfg_def.t0 = 0; % should specify this for mode 2, 3 to indicate start of
%   window
% cfg.def.t1 = 0; % should specify for mode 2, 3 to indicate end of window
%
% MvdM 2015-02-09 initial version

cfg_def.mode = 1; % 1: shuffle labels, 2: shuffle ISIs

cfg = ProcessConfig2(cfg_def,cfg_in);

switch cfg.mode
    
    case 1 % shuffle labels
   
        r = randperm(length(ts_in.t));
        ts_in.t = ts_in.t(r);
        
    case 2 % shuffle ISIs
        
        ts_in = restrict2(ts_in,cfg.t0,cfg.t1);
        
        for iC = 1:length(ts_in.t)
            
           isi = diff(cat(1,cfg.t0,ts_in.t{iC}));
           r = randperm(length(isi));
           ts_in.t{iC} = cfg.t0+cumsum(isi(r));
            
        end     
        
    case 3 % random uniform
        
        ts_in = restrict2(ts_in,cfg.t0,cfg.t1);
        
        for iC = 1:length(ts_in.t)
           
            N = length(ts_in.t{iC});
            
            ts_in.t{iC} = cfg.t0 + (cfg.t1-cfg.t0)*rand(N,1);
            
        end
end