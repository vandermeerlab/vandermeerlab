function [evt,SWR,MUA] = detector(cfg_in,SWRfreqs,csc,S)
%DETECTOR Detect candidate replay events in S and csc.
%   evt = detector(cfg_in,SWRfreqs,csc,S)
%
%   This function is intended for final analysis only. To determine input
%   parameters, work with the following functions individually in a
%   script: 
%   
%   SWRfreak()
%   amSWR()
%   amMUA()
%   precand()
%
%   Visualize output scores using MultiRaster() and sidekick().
%
%   CONFIGS (with defaults):
%   cfg.fs = 2000; in Hz, the sampling frequency
%   cfg.spkcap = 2; see amMUA
%   cfg.noisefloor = 4; see amMUA
%   cfg.mindur = 0.02; exclude events shorter than this many s
%   cfg.threshold = 8; unitless, where to cut off the detection curve
%   cfg.minCells = [];  minimum number of active cells for an event to be
%                       kept
%
%   A. Carey, Feb 2015 


%% Tell me what you're doing
tic
cfg_in.verbose = 0; % i don't want to see what the other functions have to say

cprintf(-[0 0 1],'detector: Looking for candidate replay events...');
disp(' ');

%% SWR score

disp('Working on SWR detection...')

[SWR,~,~] = amSWR(cfg_in,SWRfreqs,csc);

%% MUA score

disp('Working on MUA detection...')

[MUA,~,~] = amMUA(cfg_in,S,csc.tvec);

%% threshold

disp('Compiling and thresholding...')

evt = precand(cfg_in,csc.tvec,SWR,MUA,S);

toc

end

