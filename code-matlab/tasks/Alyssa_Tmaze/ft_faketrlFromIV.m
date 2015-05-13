function trl = ft_faketrlFromIV(cfg_in,iv_in,tvec,hdr)
%FT_FAKETRLFROMIV define 'fake' event segments of data that will be used for
% comparing results to real events (see ft_maketrlFromIV). The output is 
% deterministic and consists of a set of FieldTrip trials contained within the
% boundaries of iv_in. Note that the fake trials may overlap with one
% another in some cases (like if iv_in is short and you want a lot of fake
% trials).
%
% trl = ft_faketrlFromIV(cfg_in,iv_in,tvec,hdr)
%
% INPUTS
%
%  cfg options
%    cfg.twin: [pre post] times in s to define start and end, default [-1 1]
%    cfg.num:  the number of fake trials you want to generate (default 100)
%    cfg.mode: 'nlx' for neuralynx-based input times (default) or 'ft' for
%             fieldtrip-based (untested)
% 
%  Required inputs
%    iv_in: IV data 
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
% For more information, open ft_faketrlFromIV and read the section "about
% this function"
%
% ACarey, March 2015 TESTED FOR NLX MODE ONLY
% modified from MvdM ft_maketrl 2014-07-03

%% ABOUT THIS FUNCTION

% PREMISE
% To make sure that the fake event times do not straddle the borders of the 
% intervals, half the difference of the time window is added to the tstarts 
% and removed from the tends (see code). Then, tvec is restricted to the regions 
% contained within these intervals. Fake event center times are then evenly 
% spaced within the restricted time vector. Upon generating the output, the
% time window flanks are added to the center times, which are then
% converted to indices into the original tvec.

% WHY WAS THIS FUNCTION WRITTEN?
% if you want to compare the theta power in the LFP for the intervals of
% your detected events (SWR/MUA) to regions of the LFP likely to contain
% high theta power (like when the rat is running on the track), you can use
% this function to generate "fake events" during running times if you input
% the time intervals when the rat was on the track. The fake event trials and
% real event trials can be used with other FieldTrip functions to generate
% PSDs.
% The input iv must have been made using the iv constructor function.

% EXAMPLE USAGE for T-maze data:

%fc = ExpKeys.goodSWR(1);
%data = ft_read_neuralynx_interp(fc);
%LoadMetadata % to get trial intervals
%cfg = []; 
%cfg.num = 1000; % generate 1000 fake events (but should maybe set this equal to the number of detected real events you have?)
%cfg.twin = [-1 1]; % the time window will begin one second before and end one second after the center time of the fake event
%trl_running = ft_faketrlFromIV(cfg,metadata.taskvars.trial_iv,csc.tvec,data.hdr); 


%% parse cfg paramters
cfg_def.twin = [-1 1]; 
cfg_def.num = 100;
cfg_def.mode = 'nlx';
cfg = ProcessConfig2(cfg_def,cfg_in);

%% restrict tvec and generate tcent

% make sure the output trials are completely contained within the input iv.
% (can be done by removing the cfg.twin flanks from the intervals, because
% easier?) Also, take the indices of the interval times in tvec. 
start_idx = nearest_idx3(iv_in.tstart - cfg.twin(1),tvec);
end_idx = nearest_idx3(iv_in.tend - cfg.twin(2),tvec);

% restrict tvec using the indices

%tvecR = [];
%for iHateLoops = 1:length(start_idx)
%    keep = tvec(start_idx(iHateLoops):end_idx(iHateLoops));
%    tvecR = [tvecR; keep];
%end

tvecR = tvec(cell2mat(arrayfun(@colon,start_idx,end_idx,'UniformOutput',0)'));

% generate tcent 
stepsize = floor(length(tvecR)/(cfg.num-1)); % cfg.num - 1 because otherwise there is an extra tcent?
tcent = tvecR(1:stepsize:length(tvecR)); 

%% make the output 
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
