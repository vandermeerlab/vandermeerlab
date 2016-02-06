function lfp_tsd = LFPpower(cfg_in,lfp_tsd)
% function lfp_tsd = LFPpower(cfg,lfp_tsd)
%
% INPUT:
% tsd with filtered LFP data
%
% OUTPUT:
% tsd with power envelope obtained by Hilbert transform
%
% CFG OPTIONS with defaults:
% cfg_def.output = 'envelope'; % 'power'
% cfg.smooth = [smoothing window size] % used to smooth the power % added
% by EC 2015-05-01
%
% MvdM 2014-06-25

cfg_def = [];
cfg_def.output = 'envelope';
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % this takes fields from cfg_in and puts them into cfg

% do some checks on the data
if ~CheckTSD(lfp_tsd)
   error('LFPpower.m: CheckTSD failed.'); 
end

nSignals = size(lfp_tsd.data,1);

% process signals
for iS = 1:nSignals
    
    if cfg.verbose; fprintf('LFPpower.m: Hilberting signal %d/%d...\n',iS,nSignals); end
    
    % check for NaNs in the data; if there are, issue a warning and replace by
    % zeros
    temp_sig = lfp_tsd.data(iS,:);
    
    nan_idx = find(isnan(temp_sig));
    
    if ~isempty(nan_idx)
        fprintf('WARNING: LFPpower.m: signal %d contains NaNs (%d).\n',iS,length(nan_idx));
        temp_sig(nan_idx) = 0;
    end
    
    % obtain power
    temp_sig = hilbert(temp_sig);
    
    if isfield(cfg_in, 'smooth') 
        temp_sig = smooth(temp_sig, cfg_in.smooth);
    end
    
    switch cfg.output
        case 'power'
            temp_sig = abs(temp_sig).^2;
        case 'envelope'
            temp_sig = abs(temp_sig);
    end
    
    % reinstate NaNs and put signal back into tsd
    temp_sig(nan_idx) = NaN;
    lfp_tsd.data(iS,:) = temp_sig;
end

    

% housekeeping
lfp_tsd.cfg.history.mfun = cat(1,lfp_tsd.cfg.history.mfun,mfun);
lfp_tsd.cfg.history.cfg = cat(1,lfp_tsd.cfg.history.cfg,{cfg});