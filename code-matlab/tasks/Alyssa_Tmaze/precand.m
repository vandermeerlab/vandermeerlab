function evt = precand(cfg_in,tvec,SWR,MUA,S)
%PRECAND Compile and threshold event detection scores
%   evt = precand(cfg_in,varargin)
%     Input scores are rescaled according to their means, then combined and
%     thresholded to return event intervals.
%
%  cfg.minCells
%
% ACarey, Jan 2015.

%% Parse cfg parameters

cfg_def.mindur = 0.02; % exclude events shorter than this many s
cfg_def.threshold = 8;
cfg_def.minCells = []; % minimum number of active cells
cfg_def.verbose = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);

%% 
% in the future it might be found that certain formulas should be applied
% to the scores, like 2a + b^2 -4c, where a b c are fixed-order inputs; in
% this case should add switch cases for the formula chosen

%% talk to me

if cfg.verbose
tic
cprintf(-[0 0 1],'precand: generating pool of candidates...');
disp(' ');
end

%% compile score vector

score = sqrt(SWR.data .* MUA.data);

score = rescmean(score,0.5);

%% Threshold and get start/end times

aboveThreshold = score>cfg.threshold;
 thresh_crossings = diff(aboveThreshold); 
 start_idx = thresh_crossings == 1; 
 stop_idx = thresh_crossings == -1; 
 
 %tvec = csc.tvec; 
 tstart = tvec(start_idx);
 tend = tvec(stop_idx);
 
 %% Exclude short-duration events, and events with too few active cells
 %mindur = cfg.mindur;
 ind = zeros(size(tstart)); % this will keep track of the indices to remove
 for iInterval = length(tstart):-1:1 
     if tend(iInterval)-tstart(iInterval) < cfg.mindur % if the detected event is too short
         ind(iInterval) = iInterval;
     else
         ind(iInterval) = [];
     end
 end
 
 tstart(ind,:) = []; % remove these rows  
 tend(ind,:) = []; 
 
 if ~isempty(cfg.minCells)
     iv_in = iv(tstart,tend);
     activeCellsIV = AddNActiveCellsIV([],iv_in,S);
     nActiveCells = activeCellsIV.usr.data;
     exclude = arrayfun(@(x) x < cfg.minCells,nActiveCells);
     tstart(exclude) = [];
     tend(exclude) = [];
     nActiveCells(exclude) = [];
 end
 
 % tell user how many you found:
 message = [num2str(length(tstart)),' precandidates found.'];
 disp(message)
 
 %% Output 
evt = iv(tstart,tend);
evt.data = score;
evt.tvec = tvec;
if exist('nActiveCells','var')
    evt.nActiveCells = nActiveCells;
end

% keep a record 
evt.cfg.history.mfun = cat(1,evt.cfg.history.mfun,mfilename);
evt.cfg.history.cfg = cat(1,evt.cfg.history.cfg,{cfg});

evt.parameters = struct('threshold',cfg.threshold,'mindur',cfg.mindur,'minCells',cfg.minCells,'amSWR',SWR.parameters,'amMUA',MUA.parameters);

if cfg.verbose
toc
end
end

