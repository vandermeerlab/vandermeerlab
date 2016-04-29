function keep_iv = GetKeepIVs(cfg_in)
% function keep_iv = GetKeepIVs(cfg)
%
% returns intervals that can be kept for MotivationalT data sessions
% NOTE: if empty, means keep all
%
% for R042, excludes intervals which were not spike sorted (due to
% referencing to ground while eating)
% for R044, excludes intervals when headstage detached
%
% assumes is run from data folder
%
% MvdM 2016-04-29

cfg_def = [];
cfg_def.R044_mode = 'activity'; % 'activity', 'detachHS', 'both'
cfg_def.S = []; % can pass this to skip loading it again

cfg = ProcessConfig2(cfg_def,cfg_in);


%
[~,fn,~] = fileparts(pwd);

switch fn(1:4)
    case 'R042'
        % deal with non-spike sorted chewing intervals
        tf = FindFile('*times.mat');
        if isempty(tf)
            error('No *times.mat file found for R042');
        else
            load(tf);
            keep_iv = iv(t_start*10^-6,t_end*10^-6); % times are in neuralynx timestamps, so convert to s
        end
        
    case 'R044'
        switch cfg.R044_mode
            case 'activity'
                % deal with HS detachments -- remove extended times with no spikes
                if isempty(cfg.S)
                    please = []; please.load_questionable_cells = 1;
                    cfg.S = LoadSpikes(please);
                end
                
                cfg_Q = []; cfg_Q.dt = 1;
                Q = MakeQfromS(cfg_Q,cfg.S);
                spk_count = tsd(Q.tvec,sum(Q.data));
                
                cfg_det = []; cfg_det.threshold = 0.5; cfg_det.dcn = '>'; cfg_det.method = 'raw';
                keep_iv = TSDtoIV(cfg_det,spk_count);
           
            case 'detachHS'
                LoadMetadata;
                keep_iv = metadata.detachIV;
                
            case 'both'
                error('Not yet implemented');
                
        end % of different R044 modes
        
    case 'R050'
        keep_iv = [];
        
    case 'R064'
        keep_iv = [];
        
    otherwise
        error('Unknown rat, wrong folder?');
        
end