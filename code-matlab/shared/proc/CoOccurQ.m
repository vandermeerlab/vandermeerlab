function [out,tt_mask] = CoOccurQ(cfg_in,Q)
% function [out,tt_mask] = CoOccurQ(cfg_in,Q)
%
% obtain expected and observed co-occurrence probabilities from Q-matrix
%
% INPUTS:
%
% Q: Q-matrix [nCells x nTimeBins]
%
% OUTPUTS:
%
% out.p0: probability (fraction of time bins) each cell is active
% out.p1: expected co-occurrence under independence assumption (p(x,y) =
%  p(x)* p(y)
% out.p2: observed conditional probability p(x|y)
% out.p3: observed co-occurrence (joint) probability p(x,y)
% out.p4: z-score of p3 against shuffled data (if requested, see config below)
%
% p1 through p4 are [nCells x nCells] matrices (unless cfg.outputFormat
% specifies otherwise)
%
% tt_mask: [nCells x nCells] containing NaNs (or whatever value is specified in
% cfg.tt_repl) where cell pairs were recorded on the same tetrode
% 
% CONFIGS:
%
% cfg_def.nShuffle = 0; % 0 means no shuffle, >0 specifies number of
%  shuffles (time bins within each cell independently)
% cfg_def.num_tt = []; % if specified, specialtreat pairs from same tt (see
%   below)
% cfg_def.tt_repl = NaN; % multiply elements from same tt with thi
%
% cfg_def.outputFormat = 'matrix'; % applies to p1 through p4
%              'matrix' outputs nCells x nCells matrix
%              'vectorU' outputs the upper triangular part of the matrix as 
%                        a vector
%
% MvdM Feb 2015
% aacarey Oct 2015, added cfg.outputFormat 

cfg_def = [];
cfg_def.nShuffle = 0;
cfg_def.num_tt = []; % if specified, specialtreat pairs from same tt (see below)
cfg_def.tt_repl = NaN; % multiply elements from same tt with this
cfg_def.outputFormat = 'matrix';

cfg = ProcessConfig2(cfg_def,cfg_in);

% idea: compute active probability for each cell indepdendently
% then, for each pair, compute co-occurrence probability and test if
% different from chance (binomial test)
Q = Q.data;
nCells = size(Q,1);

% first exclude counts > 0
Q(Q > 1) = 1;

% get expected co-occurrence under independence assumption
p0 = nanmean(Q,2); % fraction of bins each cell participates in individually

% expected co-occurrence, multiply single cell probabilities
%p1 = p0*p0';
p1 = nan(nCells);
for iI = 1:nCells
    for iJ = 1:nCells
   
        p1(iI,iJ) = p0(iI).*p0(iJ);
        
    end
end
p1(logical(eye(size(p1)))) = NaN; % probability of cell co-occurring with itself is not defined


% get observed conditional probabilities (probability of cell y active
% given that cell x is active
p2 = nan(nCells);
for iI = 1:nCells
    
   % select only those rows where cell iI is active
   keep_idx = (Q(iI,:) == 1);
   Q_temp = Q(:,keep_idx);
   
   p2(iI,:) = nanmean(Q_temp,2);
   
end

% get observed co-occurrences (should be a faster way to do)
p3 = nan(nCells);
for iI = 1:nCells
    for iJ = 1:nCells
   
        p3(iI,iJ) = nanmean(Q(iI,:).*Q(iJ,:));
        
    end
end
p3(logical(eye(size(p3)))) = NaN; % probability of cell co-occurring with itself is not defined

out.p0 = p0;
out.p1 = p1;
out.p2 = p2;
out.p3 = p3;

% get co-occurrences for shuffled Q-matrix
if cfg.nShuffle > 0
    
    shuf_p4 = nan(cfg.nShuffle,nCells,nCells);
    nBins = size(Q,2);
    
    for iShuf = 1:cfg.nShuffle
        
        % permute Q-matrix
        this_Q = Q;
        for iCell = 1:nCells
            this_Q(iCell,:) = this_Q(iCell,randperm(nBins));
        end
        
        % now compute shuffled co-occurrences
        for iI = 1:nCells
            for iJ = 1:nCells
                
                shuf_p4(iShuf,iI,iJ) = nanmean(this_Q(iI,:).*this_Q(iJ,:));
                
            end
        end % of loop over cell pairs

    end % of loop over shuffles
    
    % z-score
    p4 = nan(nCells);
    
    for iI = 1:nCells
        for iJ = 1:nCells
            
            p4(iI,iJ) = (p3(iI,iJ)-nanmean(sq(shuf_p4(:,iI,iJ))))./nanstd(sq(shuf_p4(:,iI,iJ)));
            
        end
    end % of loop over cell pairs
    
    out.p4 = p4;
    
end % of shuffle conditional


%
tt_mask = [];
if ~isempty(cfg.num_tt)
   % construct same-tetrode mask
   if length(cfg.num_tt) ~= nCells
       error('tt_num has unexpected number of cells.');
   end
    
   tt_mask = ones(nCells);
   
   for iC = 1:nCells
       
      tt_mask(iC,cfg.num_tt == cfg.num_tt(iC)) = cfg.tt_repl;
       
   end
   
   out.p1 = out.p1.*tt_mask;
   out.p2 = out.p2.*tt_mask;
   out.p3 = out.p3.*tt_mask;
   
   if cfg.nShuffle > 0
       out.p4 = out.p4.*tt_mask;
   end

end

% set output format
switch cfg.outputFormat
    case 'matrix'
        % do nothing
    case 'vectorU'
        outfields = fieldnames(out);
        for iField = 2:length(outfields) % starting at 2 because we want to ignore p0, which is not a matrix
            keep = triu(ones(size(out.(outfields{iField}))),1); % get indices of values we want to keep (the 1 means we're discarding the diagonal)
            out.(outfields{iField}) = out.(outfields{iField})(find(keep));
        end
    otherwise
        error('Unrecognized cfg.outputFormat')
end
end

