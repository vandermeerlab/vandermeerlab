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
% cfg_def.InteractiveMode = 0; % if 1, use interactive mode (user input y/n
%  to swap)
%
% MvdM 2015-01-29 initial version

cfg_def.PlotOutput = 0;
cfg_def.reorder = 1; % if > 0 reversed lags, flip them
cfg_def.InteractiveMode = 0;

cfg = ProcessConfig2(cfg_def,cfg_in);

nCells = length(S_in.t);
pk_lag = nan(nCells); % make full matrix in case we want to do non-consecutive flips later
asy = nan(nCells); % make full matrix in case we want to do non-consecutive flips later

prm = 1:nCells; % new index for swapped order
swap_count = 0;

%% check place cell ordering
x = []; y = [];
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
        fh = figure;
        plot(tvec,cc,'k','LineWidth',2);
        hold on;
        if pk_lag(c1,c2) >= 0
            plot(pk_lag(c1,c2),pk_val,'.r','MarkerSize',40);
        else
            plot(pk_lag(c1,c2),pk_val,'.g','MarkerSize',40);
        end
        title(sprintf('%d-%d, lag %.2f, asy %.2f',c1,c2,pk_lag(c1,c2),asy(c1,c2)));
        
        if cfg.InteractiveMode
           
            % set TC color
           if isfield(cfg,'tch')
              set([cfg.tch(c1) cfg.tch(c2)],'FaceColor',[1 0 0]); 
           end
            
           swp = input(sprintf('%d-%d swap (y/n): ',c1,c2),'s');
           
           switch swp
               case 'y'
                   x(swap_count + 1) = c1;
                   y(swap_count + 1) = c2;
                   swap_count = swap_count + 1;
               otherwise
                   disp('Not swapped');
           end
           
           % unset TC color
           if isfield(cfg,'tch')
              set([cfg.tch(c1) cfg.tch(c2)],'FaceColor',[0 0 0]); 
           end
           
           close(fh);
           
        end
        
    end
    
end

%%
% try to improve

if ~cfg.InteractiveMode
    [x,y] = ind2sub(size(pk_lag),find(pk_lag > 0 & asy > 0));
end

vals.N = length(x);
fprintf('** There are %d consecutive pairs with positive lags & asymmetric CCFs.\n',vals.N);

S_out = S_in;

if (vals.N > 0) & cfg.reorder
    
    for iP = 1:length(x)
        
        % swap
        p1 = prm(x(iP)); p2 = prm(y(iP));
        prm(x(iP)) = p2; prm(y(iP)) = p1;
        
    end
    
    S_out.t = S_out.t(prm);
    S_out.label = S_out.label(prm);
    
end

vals.pk_lag = pk_lag;
vals.asy = asy;
vals.perm_idx = prm;