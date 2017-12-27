function boundaries = Exactness(cfg,IVann,IVdet)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


cfg_def.showFig = 1;
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

if cfg.verbose
   disp([mfun,': comparing interval boundaries...'])
end

% Restrict the IVs to manual identification times:
cfg_temp = []; cfg_temp.straddle = 1; cfg_temp.verbose = 0;
IVann = RestrictIV(cfg_temp,IVann,IVann.hdr.segments);
IVdet = RestrictIV(cfg_temp,IVdet,IVann.hdr.segments);

ratings = [1,2,3,4,5];

for iRating = 1:5
    cfg_temp = [];
    cfg_temp.verbose = 0;
    cfg_temp.threshold = iRating;
    cfg_temp.operation = '=';
    [IV,idx_0] = SelectIV(cfg_temp,IVann,'annotation');
    
    cfg_temp = []; cfg_temp.verbose = 0;
    [IV_kept,idx] = OverlapIV(cfg_temp,IV,IVdet);
    [IV_det_kept,idx0] = OverlapIV(cfg_temp,IVdet,IV);
    HitRate(iRating) = length(unique(idx))/length(unique(idx_0));

    assert(length(idx) == length(idx0),'There appears to be something wrong with OverlapIV?')
    
    % some measures for assessing "closeness" might be:
    
    % how close are the centers? (should account for splitting/merging??)
    C0 = IVcenters(IV);
    C = IVcenters(IVdet);
    
    CenterDiff = NaN(size(idx));
    
    for iID = 1:length(idx)
        CenterDiff(iID) = C0(idx(iID)) - C(idx0(iID));
    end
    boundaries.CenterDiff{iRating} = CenterDiff;
    
    % how close are the starts and ends? (should account for splitting/merging??)
    
    StartDiff = NaN(size(idx)); EndDiff = StartDiff;
    for iID = 1:length(idx)
        StartDiff(iID) = IV.tstart(idx(iID)) - IVdet.tstart(idx0(iID));
        EndDiff(iID) = IV.tend(idx(iID)) - IVdet.tend(idx0(iID));
    end
    boundaries.StartDiff{iRating} = StartDiff;
    boundaries.EndDiff{iRating} = EndDiff;
    
end

%%
if cfg.showFig
    FontSize = 14;
    
    col = linspecer(6); %{'r' 'y' 'g' 'c' 'b' 'k'};
    

    figure; set(gcf,'Name','Boundary Similarity')
    for iRating = 1:5
        
        hold on;
        scatter3(boundaries.StartDiff{1,iRating},boundaries.EndDiff{1,iRating},boundaries.CenterDiff{1,iRating},'filled','MarkerFaceColor',col(iRating,:),'MarkerEdgeColor','none')
    end
    title('Annotated minus Detected')
    zlabel('Center time differences'); xlabel('Start time differences'); ylabel('End time differences')
end

end

