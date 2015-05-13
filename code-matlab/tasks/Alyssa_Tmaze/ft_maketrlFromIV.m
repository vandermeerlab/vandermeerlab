function trl = ft_maketrlFromIV(cfg_in,iv_in,tvec,hdr)
%FT_MAKETRLFROMIV define the segments of data that will be used for
% further processing and analysis
% This is a vandermeerlab function used in place of FieldTrip's
% ft_definetrial. See ft_definetrial for additional information.
%
% trl = ft_maketrlFromIV(cfg_in,iv_in,tvec,hdr)
%
% INPUTS
%
%  cfg options
%    cfg.twin: [pre post] times in s to define start and end, e.g. [-2 2] in
%             seconds; default just uses the event time window
%    cfg.mode: 'nlx' for neuralynx-based input times (default) or 'ft' for
%             fieldtrip-based (untested)
% 
%  Required inputs
%    iv_in: IV data (like replay candidate intervals);
%     tvec: time vector (csc.tvec from LoadCSC)
%      hdr: header of ft data structure you want indices for (ex: data.hdr from
%           ft_read_neuralynx_interp)
%
% OUTPUTS
% 
%    trl:  fieldtrip-digestible [N x 3] trl variable with sample idx for start, end,
%          and trigger offset relative to start
%
%  * note that extra trial info (e.g. type) can be specified in
%   additional columns
%
% ACarey, March 2015 TESTED FOR NLX MODE ONLY
% modified from MvdM ft_maketrl 2014-07-03


cfg_def.twin = []; 
cfg_def.mode = 'nlx';
cfg = ProcessConfig2(cfg_def,cfg_in);

tcent = (iv_in.tstart+iv_in.tend)/2;
if isempty(cfg.twin)
    % then the user did not specify a different time window to use
    cfg.twin = [iv_in.tstart-tcent iv_in.tend-tcent];
end

% construct tvec based on mode
if strcmp(cfg.mode,'nlx')
    % times are on Neuralynx timebase, subtract time of first sample
 
    % convert to double; also convert FirstTimeStamp to have same units
    % as csc.tvec
    convFact = hdr.FirstTimeStamp/tvec(1);
    hdr.FirstTimeStamp = double(hdr.FirstTimeStamp/convFact);
    
    tcent = tcent - hdr.FirstTimeStamp;
    tvec = tvec-hdr.FirstTimeStamp; %minus first time to align to zero
    
else %equivalent to elseif cfg.mode = 'ft'
    tvec = cat(2,0,cumsum(repmat(1./hdr.Fs,[1 hdr.nSamples-1])));
end

% fill in trl
trl(:,3) = nearest_idx3(tcent,tvec);
trl(:,2) = nearest_idx3(tcent+cfg.twin(:,2),tvec);
trl(:,1) = nearest_idx3(tcent+cfg.twin(:,1),tvec);

trl(:,3) = trl(:,1) - trl(:,3);
