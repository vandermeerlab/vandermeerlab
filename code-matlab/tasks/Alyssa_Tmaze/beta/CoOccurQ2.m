function [out,tt_mask] = CoOccurQ2(cfg_in,Q)
% COOCCURQ2 Obtain expected and observed co-occurrence probabilities from Q-matrix 
% 
%  [cc,tt_mask] = COOCCURQ2(cfg,Q) computes the coactivation probabilities 
%    of cell pairs within the timebins in Q. Co-occurrence is a simple way to 
%    measure coordinated neural activity that does not require large neural
%    ensembles as compared to sequence analysis. Note that CoOccurQ does
%    not consider the number of times a cell spikes in a given time bin:
%    rather, it converts the Q-matrix to binary and considers whether the
%    cell was active or not, regardless of activity level.
%    
%    INPUTS:
%
%       Q: Q-matrix [nCells x nTimeBins] containing counts of the number of 
%             spikes each cell released in each time bin.
%
%    OUTPUTS:
%
%      cc: 1 x 1 struct containing different measures of co-occurrence:
%
%        .p0: [nCells x 1] double with individual values describing the 
%             probability (fraction of time bins) each cell is active.
%        .p1: [nCells x nCells] double (symmetrical matrix)* containing the 
%             expected co-occurrence under independence assumption 
%             (p(x,y) = p(x)* p(y). In other words, given cell A's
%             activity, what is the probability that we see cell B active
%             in the same time bin?
%        .p2: [nCells x nCells] double (symmetrical matrix)* containing 
%             observed conditional probability p(x|y). In other words, in
%             the bins that A fired, what proportion did B also fire?
%        .p3: [nCells x nCells] double (symmetrical matrix)* containing 
%             observed co-occurrence (joint) probability p(x,y). This is
%             the actual frequency of cell A and cell B firing together.
%        .p4: [nCells x nCells] double (symmetrical matrix)* containing 
%             the z-score of p3 against shuffled data (if requested, see 
%             config below). How much are A and B active together (in SDs
%             above the mean) as compared to randomly shuffled data. Note
%             that the number of shuffles may greatly affect the output (at
%             least 10000 shuffles are recommended for permutation tests).
%        .p5  [nCells x nCells] double (symmetrical matrix)* containing
%             the shuffled data that p4 is based on.
%
%     tt_mask: [nCells x nCells]* containing NaNs (or whatever value is 
%              specified in cfg.tt_repl) where cell pairs were recorded on
%              the same tetrode. The mask is a precaution against the
%              limitations of spike sorting.
%
%     * p1 through p4 and tt_mask are [nCells x nCells] matrices unless 
%           cfg.outputFormat specifies otherwise.
% 
%    CONFIGS:
%
%        cfg.nShuffle = 0; % 0 means no shuffle (p4 will not be computed), 
%              > 0 specifies number of shuffles (time bins within each cell 
%              independently)
%        cfg.shuffleQ = []; The Q matrix that will be used for the
%              permutation for p4. (For example, you might make shuffleQ from
%              a whole set of detected intervals, but the input Q from a
%              subset of the detected intervals).
%        cfg.useMask = 0; % If 1, disregard pairs from the same tetrode, if
%              0 don't. This requires that Q has the usr field tt_num
%              inherited from S. (Note: if tt_mask is requested as output,
%              this is automatically set to 1).
%        cfg.tt_repl = NaN; % Replace elements from same tt with this. This
%              config field applies only if cfg.useMask = 1.
%        cfg.outputFormat = 'matrix'; % Applies to p1 through p4
%              'matrix' outputs nCells x nCells matrix
%              'vectorU' outputs the upper triangular part of the matrix as 
%                        a vector
%        cfg.resetRNG = 0; If 1, temporarily resets the rng settings so that
%               different calls will produce the same output for p4 (the
%               permutation test). After running, the original rng settings
%               are reinstated. If 0, the current rng state is used and
%               repeated runs with the same cfg and inputs will not produce
%               the same output. This config option is useful for testing.
%
%  See also: MakeQfromS2
%
% MvdM Feb 2015, initial version
% aacarey Dec 2015, version 2

cfg_def = [];
cfg_def.resetRNG = 0; % if 1, temporarily resets rng
cfg_def.nShuffle = 0;
cfg_def.shuffleQ = [];
cfg_def.useMask = 0; % if specified, specialtreat pairs from same tt (see below)
cfg_def.tt_repl = NaN; % multiply elements from same tt with this
cfg_def.outputFormat = 'matrix';
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);


if cfg.verbose; fprintf('\n%s: computing coactivation probabilities...\n',mfun); end

if nargout > 1; cfg.useMask = 1; end

if cfg.resetRNG; originalRNGsettings = rng; rng(1,'twister'); end

% idea: compute active probability for each cell indepdendently
% then, for each pair, compute co-occurrence probability and test if
% different from chance (binomial test)

if cfg.useMask; tt_num = Q.usr.tt_num; end
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
    if cfg.verbose; disp([mfun,': shuffling Q-matrix ',num2str(cfg.nShuffle),' times...']); end
    if isempty(cfg.shuffleQ)
        Qshuf = Q;
    else
        Qshuf = cfg.shuffleQ.data;
        Qshuf(Qshuf > 1) = 1;
    end
    nRows = size(Qshuf,1);
    
    shuf_p4 = nan(cfg.nShuffle,nRows,nRows);
    nBins = size(Qshuf,2);
    
    for iShuf = 1:cfg.nShuffle
        
        % permute Q-matrix
        this_Q = Qshuf;
        for iCell = 1:nRows
            this_Q(iCell,:) = this_Q(iCell,randperm(nBins));
        end
        
        % now compute shuffled co-occurrences
        for iI = 1:nRows
            for iJ = 1:nRows
                
                shuf_p4(iShuf,iI,iJ) = nanmean(this_Q(iI,:).*this_Q(iJ,:));
                
            end
        end % of loop over cell pairs

    end % of loop over shuffles
    
    % z-score
    p4 = nan(nRows);
    
    for iI = 1:nRows
        for iJ = 1:nRows
            
            p4(iI,iJ) = (p3(iI,iJ)-nanmean(sq(shuf_p4(:,iI,iJ))))./nanstd(sq(shuf_p4(:,iI,iJ)));
            
        end
    end % of loop over cell pairs
    
    out.p4 = p4;
    out.p5 = shuf_p4;
    
end % of shuffle conditional


% handle tt_mask (for cells that were recorded on the same tt)
if cfg.useMask
    % construct same-tetrode mask
    if length(tt_num) ~= nCells
        error('Q is poorly formed: Q.usr.tt_num must have the same number of entries as Q.data has rows.');
    end
    
    tt_mask = ones(nCells);
    
    for iC = 1:nCells
        
        tt_mask(iC,tt_num == tt_num(iC)) = cfg.tt_repl;
        
    end
    
    out.p1 = out.p1.*tt_mask;
    out.p2 = out.p2.*tt_mask;
    out.p3 = out.p3.*tt_mask;
    
    if cfg.nShuffle > 0
        out.p4 = out.p4.*tt_mask;
    end
end % of tt_num empty check


% set output format
switch cfg.outputFormat
    case 'matrix'
        % do nothing
    case 'vectorU'
        outfields = fieldnames(out);
        keep = triu(ones(nCells),1); % get indices of values we want to keep (the 1 means we're discarding the diagonal)
        for iField = 2:length(outfields) % starting at 2 because we want to ignore p0, which is not a matrix
            out.(outfields{iField}) = out.(outfields{iField})(find(keep));
        end
        if nargout > 1
            tt_mask = tt_mask(find(keep)); % also do this for tt_mask so it matches the output
        end
    otherwise
        error('Unrecognized cfg.outputFormat option specified')
end

if cfg.resetRNG; rng(originalRNGsettings); end
end

