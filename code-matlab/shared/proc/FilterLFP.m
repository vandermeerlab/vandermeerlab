function lfp_tsd = FilterLFP(cfg_in,lfp_tsd)
% function lfp_tsd = FilterLFP(cfg,lfp_tsd)
%
% INPUT:
% tsd with LFP data
%
% OUTPUT:
% tsd with filtered LFP data
%
% CFG OPTIONS with defaults:
%
% cfg.type = 'butter'; % {'cheby1','butter','fdesign'} -- type of filter to be used
% cfg.order = 4; % filter order
% cfg.display_filter = 0; % show output of fvtool on filter
% cfg.bandtype = 'bandpass'; % 'highpass', 'lowpass'
% cfg.R = 0.5; % passband ripple (in dB) for Chebyshev filters only
% cfg.f = [6 10]; filter range to use (in Hz)
% cfg.verbose = 1; If 1 display helpful text in command window, if 0 don't
%
% MvdM 2014-06-24, 25 (update cfg_in)

cfg_def.type = 'butter';
cfg_def.order = 4;
cfg_def.display_filter = 0;
cfg_def.band = 'bandpass'; % not used
cfg_def.R = 0.5; % passband ripple (in dB) for Chebyshev filters only
cfg_def.f = [6 10];
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % this takes fields from cfg_in and puts them into cfg

% do some checks on the data
if ~CheckTSD(lfp_tsd)
   error('FilterLFP.m: CheckTSD failed.'); 
end

% check reported Fs in headers
nSignals = length(lfp_tsd.cfg.hdr);
for iS = nSignals:-1:1
    
   if isfield(lfp_tsd.cfg.hdr{iS},'SamplingFrequency')
       reported_Fs(iS) = lfp_tsd.cfg.hdr{iS}.SamplingFrequency;
   elseif isfield(lfp_tsd.cfg.hdr{iS},'Fs')
       reported_Fs(iS) = lfp_tsd.cfg.hdr{iS}.Fs;
   else
      error('Unknown Fs.'); 
   end
   
    
end

reported_Fs = unique(reported_Fs);
if length(reported_Fs) > 1
   error('FilterLFP.m: multiple sampling frequencies in header.'); 
end

tvec_diffs = diff(lfp_tsd.tvec);
median_Fs = 1./median(tvec_diffs);

if cfg.verbose; fprintf('FilterLFP.m: reported Fs %.2f, median tvec Fs %.2f.\n',reported_Fs,median_Fs); end

Fs = reported_Fs;

% construct the filter
Wn = cfg.f ./ (Fs/2); % convert to units of half sample rate
switch cfg.type
   
    case 'butter'
        
        %[z,p,k] = butter(cfg.order,Wn);
        [b,a] = butter(cfg.order,Wn);
        
    case 'cheby1'
        
        %[z,p,k] = cheby1(cfg.order,cfg.R,Wn);
        [b,a] = cheby1(cfg.order,cfg.R,Wn);
    
    case 'fdesign'
        
        d = fdesign.bandpass('N,F3dB1,F3dB2',cfg.order,cfg.f(1),cfg.f(2),Fs);
        Hd = design(d,'butter');
        b = Hd.sosMatrix; a = Hd.scaleValues;
end

% convert to SOS format
%[sos,g] = zp2sos(z,p,k);

% display if requested
if cfg.display_filter
    %h = dfilt.df2sos(sos,g);
    %fvtool(h);
    if exist('Hd','var')
        fvtool(Hd);
    else
        fvtool(b,a);
    end
    fprintf('FilterLFP.m: paused, press key to continue...\n');
    pause;
end

% process signals
for iS = 1:nSignals
    
    if cfg.verbose; fprintf('FilterLFP.m: filtering signal %d/%d...\n',iS,nSignals); end
    
    % check for NaNs in the data; if there are, issue a warning and replace by
    % zeros
    temp_sig = lfp_tsd.data(iS,:);
    
    nan_idx = find(isnan(temp_sig));
    
    if ~isempty(nan_idx)
        fprintf('WARNING: FilterLFP.m: signal %d contains NaNs (%d).\n',iS,length(nan_idx));
        temp_sig(nan_idx) = 0;
    end
    
    % filter
    %temp_sig = filtfilt(sos,g,temp_sig);
    temp_sig = filtfilt(b,a,temp_sig);
    
    % reinstate NaNs and put signal back into tsd
    temp_sig(nan_idx) = NaN;
    lfp_tsd.data(iS,:) = temp_sig;
end

% housekeeping
lfp_tsd.cfg.history.mfun = cat(1,lfp_tsd.cfg.history.mfun,mfun);
lfp_tsd.cfg.history.cfg = cat(1,lfp_tsd.cfg.history.cfg,{cfg});