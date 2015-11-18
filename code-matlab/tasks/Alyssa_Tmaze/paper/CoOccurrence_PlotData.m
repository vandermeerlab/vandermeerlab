%% config
cfg = [];
cfg.rats = {'R042','R044','R050','R064'};

clear cooc_results
%% COLLECT CO-OCCURRENCE DATA INTO INDIVDUAL PVALS FOR EACH RAT (for plotting
% later)
COMBp.p0 = []; % combined pvals for all food sessions for all rats;
COMBp.p4 = [];

cooc_results.all.foodL = COMBp; % left cells on food day
cooc_results.all.foodR = COMBp; % right cells on food day
cooc_results.all.waterL = COMBp;
cooc_results.all.waterR = COMBp;

pnames = fieldnames(COMBp); % for the pval concatenating loop down below
p0_idx = strmatch('p0',pnames); % make sure we know which idx corresponds to p0, which needs to be handled differently when averaging

for iRat = 1:length(cfg.rats)
    disp(' '); disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['~~  COLLECTING CO-OCCURRENCE DATA FOR ',cfg.rats{iRat},' ~~'])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    % initialize combined p vals struct
    %pvalFields = fieldnames(cooc.(cfg.rats{iRat})(iSession).(arms{iArm}).ALLp); % pvalFields will be 'p0','p1','p2' and so on
    %COMBp = cell2struct(cell(size(pvalFields)),pvalFields,1); % this takes the field names, generates an empty cell, and then turns the cell into an empty struct with those field names
    
    cooc_results.(cfg.rats{iRat}).foodL = COMBp;
    cooc_results.(cfg.rats{iRat}).foodR = COMBp;
    cooc_results.(cfg.rats{iRat}).waterL = COMBp;
    cooc_results.(cfg.rats{iRat}).waterR = COMBp;
    
    for iSession = 1:length(cooc_data.(cfg.rats{iRat}))
        
        restriction = cooc_data.(cfg.rats{iRat})(iSession).restrictionType;
        for iArm = 1:length(arms)
            
            for iPname = 1:length(pnames)
                
                this_coOccur = cooc_data.(cfg.rats{iRat})(iSession).(arms{iArm}).ALLp.(pnames{iPname}); % this is a nCells x nCells symmetric matrix
                
                if iPname ~= p0_idx
                    this_coOccur_mask = triu(ones(size(this_coOccur)),1); % only take upper half of the matrix
                    this_coOccur = this_coOccur(logical(this_coOccur_mask));
                end
                
                % add these data to collection
                cooc_results.(cfg.rats{iRat}).([restriction,arms{iArm}]).(pnames{iPname}) = cat(1,cooc_results.(cfg.rats{iRat}).([restriction,arms{iArm}]).(pnames{iPname}),this_coOccur);
                cooc_results.all.([restriction,arms{iArm}]).(pnames{iPname}) = cat(1,cooc_results.all.([restriction,arms{iArm}]).(pnames{iPname}),this_coOccur);
            
            end % fields
        end % arms
    end % sessions
end % rats
