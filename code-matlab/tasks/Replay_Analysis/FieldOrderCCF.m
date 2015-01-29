function [S_out,vals] = FieldOrderCCF(cfg_in,S_in)
% function [S_out,vals] = FieldOrderCCF(cfg_in,S_in)
%
% reorder ordered place fields based on cross-correlation peak lags: the
% expectation is that for correctly ordered fields, each consecutive pair
% has a negative lag ccf peak (Skaggs et al. 1996; Dragoi et al. 2006)
%
% currently does a single pass in which pairs with reverse lags are flipped
%
% INPUTS:
%
% S_in: ts struct, pre-ordered to that S_in.t{1} is the first place cell
% (such as is the case for the output of MakeTC()) 
%
% OUTPUTS:
%
% S_out: reordered ts struct
% vals.N: count of remaining pairs for which the ccf lag is reversed
% vals.pk_lag: matrix of ccf peak lags (in s)
% vals.asy: matrix of ccf asymmetries (contrast between positive and
%  negative lag halves)
% vals.perm_idx: S_in.t(perm_idx) = S_out.t
%
% CFG:
%
% cfg.PlotOutput = 0;
% cfg_def.reorder = 1; % if > 0 reversed lags, flip them
%
% MvdM 2015-01-29 initial version

cfg_def.PlotOutput = 0;
cfg_def.reorder = 1; % if > 0 reversed lags, flip them
cfg = ProcessConfig2(cfg_def,cfg_in);

nCells = length(S_in.t);
pk_lag = nan(nCells); % make full matrix in case we want to do non-consecutive flips later
asy = nan(nCells); % make full matrix in case we want to do non-consecutive flips later

%% check place cell ordering
for iC = 1:nCells-1

    c1 = iC; c2 = c1 + 1;
    
    ts1 = S_in.t{c1}; ts2 = S_in.t{c2};
    
    cfg_ccf = [];
    cfg_ccf.smooth = 1;
    [cc,tvec] = ccf(cfg_ccf,ts1,ts2);
    
    [pk_val,pk_idx] = max(cc);
    pk_lag(c1,c2) = tvec(pk_idx);
    
    asy_l = nanmean(cc(1:floor(length(cc)/2)));
    asy_r = nanmean(cc(ceil(length(cc)/2)+1:end));
    asy(c1,c2) = (asy_r-asy_l)./(asy_l+asy_r);
    
    
    if cfg.PlotOutput
        figure;
        plot(tvec,cc,'k','LineWidth',2);
        hold on;
        if pk_lag(c1,c2) >= 0
            plot(pk_lag(c1,c2),pk_val,'.r','MarkerSize',40);
        else
            plot(pk_lag(c1,c2),pk_val,'.g','MarkerSize',40);
        end
        title(sprintf('%d-%d, lag %.2f, asy %.2f',c1,c2,pk_lag(c1,c2),asy(c1,c2)));
    end
    
end

%%
% try to improve

[x,y] = ind2sub(size(pk_lag),find(pk_lag > 0 & asy > 0));
vals.N = length(x);
vals.pk_lag = pk_lag;

vals.asy = asy;

fprintf('** There are %d consecutive pairs with positive lags & asymmetric CCFs.\n',vals.N);

S_out = S_in;
prm = 1:nCells;

if (vals.N > 0) & cfg.reorder
    
    for iP = 1:length(x)
        
        % swap
        p1 = prm(x(iP)); p2 = prm(y(iP));
        prm(x(iP)) = p2; prm(y(iP)) = p1;
        
    end
    
    S_out.t = S_out.t(prm);
    S_out.label = S_out.label(prm);
    
end

vals.perm_idx = prm;