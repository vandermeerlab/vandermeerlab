function ROCdata = GetROCdata(cfg,IVann,TSD)
%GETROCDATA Summary of this function goes here
%   Detailed explanation goes here
%
%   INPUTS
%       cfg
%       ivANN - annotated IV data, has usr.annotation containing
%                human's confidence rating
%       TSD   - tsd struct containing detector data, such as the output
%               from OldWizard or any such SWR detection algorithm
%
% aacarey NOV 2017

cfg_def.thresholds = 1:0.1:10;
cfg_def.method = 'zscore';
cfg_def.verbose = 1;

mfun = mfilename;

cfg = ProcessConfig(cfg_def,cfg,mfun);

% make sure evt is annotation
if ~isfield(IVann.usr,'annotation')
    error('IVann data is not annotated; data must be annotated')
end

% Restrict IVann to definite segment intervals
cfg_temp = []; cfg_temp.verbose = 0;
IVann = RestrictIV(cfg_temp,IVann,IVann.hdr.segments.tstart,IVann.hdr.segments.tend);

% allocate space
hitrates1 = nan(1,length(cfg.thresholds));
hitrates2 = hitrates1; hitrates3 = hitrates1; hitrates4 = hitrates1; hitrates5 = hitrates1; falseposrates = hitrates1;
nEvt1 = hitrates1; nEvt2 = hitrates1; nEvt3 = hitrates1; nEvt4 = hitrates1; nEvt5 = hitrates1; nFalsePos = hitrates1;

% make intervals from different thresholds
for iThr = length(cfg.thresholds):-1:1

    cfg_temp = []; cfg_temp.threshold = cfg.thresholds(iThr); cfg_temp.method = 'zscore'; cfg_temp.verbose = 0;
    %cfg.ResizeAtMean = 0;
    IVatThreshold(iThr)  = TSDtoIV2(cfg_temp,TSD);
    
    cfg_temp = []; cfg_temp.verbose = 0;
    IVatThreshold(iThr) = RestrictIV(cfg_temp,IVatThreshold(iThr),IVann.hdr.segments.tstart,IVann.hdr.segments.tend);
    
    cfg_temp = []; cfg_temp.showFig = 0; cfg_temp.verbose = 0;
    [hitrate,nEvt,~] = HitRate(cfg_temp,IVann,IVatThreshold(iThr));
    
    % Collect HitRate data
    if length(hitrate)>6; error('Currently cannot handle more than 6 HitRate categories'); end
    hitrates1(iThr) = hitrate(1);
    hitrates2(iThr) = hitrate(2);
    hitrates3(iThr) = hitrate(3);
    hitrates4(iThr) = hitrate(4);
    hitrates5(iThr) = hitrate(5);
    falseposrates(iThr) = hitrate(6);
    
    % Collect nEvt data
    nEvt1(iThr) = nEvt(1);
    nEvt2(iThr) = nEvt(2);
    nEvt3(iThr) = nEvt(3);
    nEvt4(iThr) = nEvt(4);
    nEvt5(iThr) = nEvt(5);
    nFalsePos(iThr) = nEvt(6);
end

% Make output
ROCdata.hitrates1 = hitrates1;
ROCdata.hitrates2 = hitrates2;
ROCdata.hitrates3 = hitrates3;
ROCdata.hitrates4 = hitrates4;
ROCdata.hitrates5 = hitrates5;
ROCdata.falseposrates = falseposrates; % false positives
ROCdata.nEvt1 = nEvt1;
ROCdata.nEvt2 = nEvt2;
ROCdata.nEvt3 = nEvt3;
ROCdata.nEvt4 = nEvt4;
ROCdata.nEvt5 = nEvt5;
ROCdata.nFalsePos = nFalsePos;

% Record config history
ROCdata = History(ROCdata,mfun,cfg);

end