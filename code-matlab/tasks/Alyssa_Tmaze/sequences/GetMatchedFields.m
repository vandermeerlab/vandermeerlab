function [TC_leftM,TC_rightM] = GetMatchedFields(cfg_in,TC_left,TC_right)
% function [TC_leftM,TC_rightM] = GetMatchedFields(cfg_in,TC_left,TC_right)
%
% returns pair of tuning curve sets (left, right) with same number of cells
%
% algorithm: each cell in smallest set gets greedily matched to the cell in
% the largest set with the closest peak firing rate
%
% NOTE: ignores multiple fields per cell -- based on peak location only
%
% MvdM 2015-02-19 initial version


cfg_def = [];

cfg = ProcessConfig2(cfg_def,cfg_in);

nLcells = length(TC_left.template_idx);
nRcells = length(TC_right.template_idx);

match = '';
if nLcells == nRcells
    fprintf('Cell numbers equal (%d), TCs returned as is.\n',nLcells);
    TC_leftM = TC_left;
    TC_rightM = TC_right;
    return;
elseif nLcells > nRcells
    fprintf('nLcells %d, nRcells %d...\n',nLcells,nRcells);
    match = 'R';
    
    min_cell_idx = TC_right.template_idx;  % idx into S.t
    min_cell_loc = TC_right.peak_idx; % location of peak (NOTE: ignores multiple fields) 
    
    max_cell_idx = TC_left.template_idx;
    max_cell_loc = TC_left.peak_idx;
    
else
    fprintf('nLcells %d, nRcells %d...\n',nLcells,nRcells);
    match = 'L';
    
    min_cell_idx = TC_left.template_idx;  % idx into S.t
    min_cell_loc = TC_left.peak_idx; % location of peak (NOTE: ignores multiple fields) 
    
    max_cell_idx = TC_right.template_idx;
    max_cell_loc = TC_right.peak_idx;
end

%
keep_idx = 1:length(max_cell_loc); % idx of cells to keep in bigger TC
for iC = 1:length(min_cell_loc)
    
    this_loc = min_cell_loc(iC);
    
    this_matched_idx = nearest_idx3(this_loc,max_cell_loc); % idx into max_cell_loc
   
    %matched_idx(iC) = max_cell_idx(this_matched_idx);
    matched_loc(iC) = max_cell_loc(this_matched_idx);
    matched_cvt(iC) = keep_idx(this_matched_idx);
    
    max_cell_loc(this_matched_idx) = []; % remove matched trial from available pool to prevent duplicates
    %max_cell_idx(this_matched_idx) = [];
    keep_idx(this_matched_idx) = [];
    
end

matched_cvt = sort(matched_cvt,'ascend');

switch match
    case 'L'
        TC_leftM = TC_left;
        
        TC_rightM = TC_right;
        TC_rightM.template_idx = TC_rightM.template_idx(matched_cvt);
        TC_rightM.peak_idx = TC_rightM.peak_idx(matched_cvt);
        TC_rightM.peak_loc = TC_rightM.peak_loc(matched_cvt);
        
        % now also do fields
        fpeak_val = []; fpeak_idx = [];
        iPeak = 1;
        for iP = 1:length(TC_rightM.peak_loc)
            fpeak_val(iPeak:iPeak+length(TC_rightM.peak_loc{iP})-1) = TC_rightM.peak_loc{iP};
            fpeak_idx(iPeak:iPeak+length(TC_rightM.peak_loc{iP})-1) = TC_rightM.template_idx(iP);
            
            iPeak = iPeak + length(TC_rightM.peak_loc{iP});
        end
        [fpeak_val,sort_idx] = sort(fpeak_val,'ascend');
        fpeak_idx = fpeak_idx(sort_idx);
        
        TC_rightM.field_template_idx = fpeak_idx;
        TC_rightM.field_loc = fpeak_val;
        
    case 'R'
        TC_rightM = TC_right;
        
        TC_leftM = TC_left;
        TC_leftM.template_idx = TC_leftM.template_idx(matched_cvt);
        TC_leftM.peak_idx = TC_leftM.peak_idx(matched_cvt);
        TC_leftM.peak_loc = TC_leftM.peak_loc(matched_cvt);
        
        % now also do fields
        fpeak_val = []; fpeak_idx = [];
        iPeak = 1;
        for iP = 1:length(TC_leftM.peak_loc)
            fpeak_val(iPeak:iPeak+length(TC_leftM.peak_loc{iP})-1) = TC_leftM.peak_loc{iP};
            fpeak_idx(iPeak:iPeak+length(TC_leftM.peak_loc{iP})-1) = TC_leftM.template_idx(iP);
            
            iPeak = iPeak + length(TC_leftM.peak_loc{iP});
        end
        [fpeak_val,sort_idx] = sort(fpeak_val,'ascend');
        fpeak_idx = fpeak_idx(sort_idx);
        
        TC_leftM.field_template_idx = fpeak_idx;
        TC_leftM.field_loc = fpeak_val;
end
