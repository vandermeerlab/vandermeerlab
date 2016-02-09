function trl = ft_maketrl(cfg_in)
% function trl = ft_maketrl(cfg)
%
% INPUTS
%
% cfg.t: trial center times (zero or trigger points) in seconds
% cfg.twin: [pre post] times in s to define start and end, e.g. [-2 2] in
% seconds
% cfg.mode: 'ft' for fieldtrip-based input times (default) or 'nlx' for
% neuralynx-based
% cfg.hdr: hdr of ft data structure you want indices for
%
%
% OUTPUTS
% 
% fieldtrip-digestible [N x 3] trl variable with sample idx for start, end,
% and trigger offset relative to start
%
% note that extra colums with trial info (e.g. type) can be specified in
% further columns
%
% MvdM 2014-07-03

cfg_def.twin = [-2 2];
cfg_def.mode = 'ft';
cfg = ProcessConfig2(cfg_def,cfg_in);

if ~isfield(cfg,'t')
    error('cfg.t (trial times) must be specified.');
end


if ~isfield(cfg,'hdr')
    error('cfg.hdr (cfg.hdr of ft data structure to get idxs for) must be specified.');
end
    
% construct tvec based on mode
if strcmp(cfg.mode,'nlx')
    % times are on Neuralynx timebase, subtract time of first sample
    cfg.t = cfg.t - double(cfg.hdr.FirstTimeStamp)/10^6;
    %tvec = (cfg.hdr.FirstTimeStamp:cfg.hdr.TimeStampPerSample:cfg.hdr.LastTimeStamp)-double(cfg.hdr.FirstTimeStamp); %minus first time to align to zero
    tvec = cat(2,0,cumsum(repmat(1./cfg.hdr.Fs,[1 cfg.hdr.nSamples-1])));
else %equivalent to elseif cfg.mode = 'ft'
    tvec = cat(2,0,cumsum(repmat(1./cfg.hdr.Fs,[1 cfg.hdr.nSamples-1])));
end

% fill in trl
trl(:,3) = nearest_idx3(cfg.t,tvec);
trl(:,2) = nearest_idx3(cfg.t+cfg.twin(2),tvec);
trl(:,1) = nearest_idx3(cfg.t+cfg.twin(1),tvec);

trl(:,3) = trl(:,1) - trl(:,3);
