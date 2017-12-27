function [hitrate,nEvt,score] = HitRate(cfg,IVann,IVdet)
%HITRATE Measure detection success in terms of hit rate
%    hitrate = HITRATE(cfg,iv_ann,iv_det) calculates the proportion of
%    manually-identified intervals that were detected automatically. The
%    calculation is performed for each group of annotations: 1's, 2's, 3's,
%    4's, and 5's. Additionally, a false positive rate is reported.
%
%    INPUTS
%          cfg: config struct with fields controlling function behavior
%        IVann: iv data with ratings contained in a usr annotation field (see
%               AnnotateIV) - manually-identified intervals
%        IVdet: detected intervals
%    
%    OUTPUTS
%    hitrate:
%       hitrate(1): proportion of 1's that were detected
%       hitrate(2): proportion of 2's that were detected
%       hitrate(3): proportion of 3's that were detected
%       hitrate(4): proportion of 4's that were detected
%       hitrate(5): proportion of 5's that were detected
%       hitrate(6): false positive rate (intervals that were present in the
%                   detected set but not in the manually-identified set)
%
%    CONFIG OPTIONS
%     cfg.verbose = 1; Print hit rate data to the command window as
%                      percentage detected
%     cfg.showFig = 1; Display a figure of hit rate data
%
% aacarey Nov 2015
% aacarey Nov 2017

cfg_def.showFig = 1;
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

if cfg.verbose
   disp('HitRate: calculating detection success')
end

% Restrict the IVs to manual identification times:
cfg_temp = []; cfg_temp.straddle = 1; cfg_temp.verbose = 0;
IVann = RestrictIV(cfg_temp,IVann,IVann.hdr.segments);
IVdet = RestrictIV(cfg_temp,IVdet,IVann.hdr.segments);

% Calculate hit rate
for iRating = 1:5
    cfg_temp = [];
    cfg_temp.verbose = 0;
    cfg_temp.threshold = iRating;
    cfg_temp.operation = '=';
    % Select all the manually-identified ripples with a given rating (iRating)
    [IVann_iRating,idx_iRating] = SelectIV(cfg_temp,IVann,'annotation');
    
    % Check whether the human-identified region is represented within the
    % detected intervals
    cfg_temp = []; cfg_temp.verbose = 0;
    [~,idx_correct] = OverlapIV(cfg_temp,IVann_iRating,IVdet); % Get the indices of human intervals that overlap with detected intervals
    
    hitrate(iRating) = length(unique(idx_correct))/length(unique(idx_iRating));
    nEvt(iRating) = numel(idx_correct);

    % Make sure something isn't terribly wrong
    [~,idx_check] = OverlapIV(cfg_temp,IVdet,IVann_iRating); % Get the indices of detected intervals that overlap with human intervals
    assert(length(idx_correct) == length(idx_check),'There appears to be something wrong with OverlapIV?')
      
end

%% Count false positives

% False positives are any intervals in the detected set that do not overlap
% with any of the intervals in the manually-identified set
cfg_temp = []; cfg_temp.verbose = 0;
iv_flsps = DifferenceIV(cfg_temp,IVdet,IVann);

hitrate(6) = length(iv_flsps.tstart)/length(IVdet.tstart);
nEvt(6) = length(iv_flsps.tstart);
% hitrate(6) = length(iv_flsps.tstart)/length(IVann.tstart);

%% Assign score
score = hitrate(1)*5 + hitrate(2)*4 + hitrate(3)*3 + hitrate(4)*2 + hitrate(5)*1 - hitrate(6)*6;

%% Plot stuff
if cfg.showFig
    FontSize = 14;
    
    location = [1 2 3 4 5 6];
    col = linspecer(6); %{'r' 'y' 'g' 'c' 'b' 'k'};
    figure; set(gcf,'Name','Hit Rate')
    for iBar = 1:length(hitrate)
        h(iBar) = bar(location(iBar),hitrate(iBar),'FaceColor',col(iBar,:),'EdgeColor','k');
        hold on
    end
    
    title(['Score: ', num2str(score)])
    ylabel('Success rate','FontSize',FontSize') 
    xlabel('Rating','FontSize',FontSize)
    set(gca,'YLim',[0 1.01],'XTick',location,'XTickLabel',{'1' '2' '3' '4' '5' 'flsps'},'FontSize',FontSize,'LineWidth',1,'YTick',0:0.1:1);
    box off
end

% Talk to me
if cfg.verbose
   disp('Rating      Hit Rate')
   fprintf('1           %.2f%%\n',hitrate(1)*100)
   fprintf('2           %.2f%%\n',hitrate(2)*100)
   fprintf('3           %.2f%%\n',hitrate(3)*100)
   fprintf('4           %.2f%%\n',hitrate(4)*100)
   fprintf('5           %.4f%%\n',hitrate(5)*100)
   fprintf('False pos:  %.4f%%\n',hitrate(6)*100)
   fprintf('score:      %.2f\n',score)
   
end

% Keep a record
% hitrate = History(hitrate,mfun,cfg);
end

