function iv_in = CountCycles(cfg_in,tsd_in,iv_in)
% function CountCycles(cfg,tsd_in,iv_in)
%
% adds .nCycles usr field to intervals iv_in, counting the number of oscillation
% cycles in each interval; also some undocumented usr fields
%
% filters the incoming signal (tsd_in), then detects the number of peaks
% within each interval.
%
% cfg_def = [];
% cfg_def.threshold_type = 'zscore'; % {'zscore','raw'}
% cfg_def.threshold = 2; % peaks need to be over this threshold to be counted
% cfg_def.debug = 0;
% cfg_def.filter_cfg = []; % use this to do the filtering
%
% MvdM 2016-01-07 initial version

cfg_def = [];
cfg_def.threshold_type = 'zscore'; % {'zscore','raw'}
cfg_def.threshold = 2; % peaks need to be over this threshold to be counted
cfg_def.debug = 0;
cfg_def.filter_cfg = []; % use this to do the filtering

cfg = ProcessConfig(cfg_def,cfg_in);

if isempty(cfg.filter_cfg)
    error('Must specify cfg_in.filter_cfg!');
end

% filter
tsd_f = FilterLFP(cfg.filter_cfg,tsd_in);

% first, hilbertize
nan_idx = find(isnan(tsd_f.data));
tsd_f.data(nan_idx) = 0;
tsd_h = angle(hilbert(tsd_f.data));

% convert to zscore if needed
switch cfg.threshold_type
    case 'zscore'
   
        %tsd_in.data(nan_idx) = NaN;
        tsd_f.data = zscore(tsd_f.data); % isn't ideal because affected by number of zeros
        
end

% loop over iv's
nIV = length(iv_in.tstart);

for iIV = nIV:-1:1
    
   this_data = tsd_f.data(TSD_getidx2(tsd_f,iv_in.tstart(iIV),iv_in.tend(iIV)));
   this_data_raw = tsd_in.data(TSD_getidx2(tsd_in,iv_in.tstart(iIV),iv_in.tend(iIV)));
   this_data_h = tsd_h(TSD_getidx2(tsd_f,iv_in.tstart(iIV),iv_in.tend(iIV)));
   
   % peaks are zero crossings in hilbert
   temp_h = this_data_h > 0;
   this_pk_idx = find(diff(temp_h) == 1);
   this_pk_val = this_data(this_pk_idx); this_pk_val_raw = this_data_raw(this_pk_idx);
   
   % troughs are pi crossings in hilbert
   this_tr_idx = find(diff(this_data_h) < 0);
   this_tr_val = this_data(this_tr_idx); this_tr_val_raw = this_data_raw(this_tr_idx);
   
   if cfg.debug
    
       subplot(211)
       plot(this_data)
       hold on;
       plot(this_data_h/1000,'r')
       %plot(this_data_h,'r')
       
       plot(this_pk_idx,this_pk_val,'.r');
       plot(this_tr_idx,this_tr_val,'.g');
       
       subplot(212)
       plot(this_data_raw)
       hold on;
       plot(this_pk_idx,this_pk_val_raw,'.r');
       plot(this_tr_idx,this_tr_val_raw,'.g');
       
   end
   
   % track number of cycles
   iv_in.usr.nCycles(iIV) = sum(this_pk_val > cfg.threshold);
   
   % track mean filtered peak amplitude
   iv_in.usr.mean_filt(iIV) = nanmean(this_pk_val);
   iv_in.usr.min_filt(iIV) = min(this_pk_val);
   
   % track variability in how far from mean peaks and throughs are
   all_peak_tr = cat(2,this_pk_val-nanmean(this_data),this_tr_val-nanmean(this_data));
   all_peak_tr = abs(all_peak_tr);
   iv_in.usr.var(iIV) = std(all_peak_tr)./nanmean(all_peak_tr);
   
   all_peak_tr = cat(2,this_pk_val_raw-nanmean(this_data_raw),this_tr_val_raw-nanmean(this_data_raw));
   all_peak_tr = abs(all_peak_tr);
   iv_in.usr.var_raw(iIV) = std(all_peak_tr)./nanmean(all_peak_tr);

       
end