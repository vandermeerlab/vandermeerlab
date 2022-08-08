function tc = TuningCurvesIG(cfg_in,S,tuning_var)
% function tc = TuningCurvesIG(cfg,S,tuning_var)
%
% computes tuning curves for cells in S relative to
% variable(s) tuning_var
%
% INPUTS:
%
% S: ts with spike data
% tuning_var: tsd with nDim x nSamples tuning variable (e.g. position)
%
% OUTPUTS:
%
% tc.tc: nCells x nBins tuning curves
% tc.occ: 1 x nBins occupancy (sample count)
% tc.tc2D: nCells x nXBins x NYBins tuning curves (for 2D tuning_var only)
% tc.binsUsed: idx into tc.tc2D specifying which bins are used to make
%   tc.tc
% tc.pos_idx = 1 x nSamples idx, assignment of tuning variable(s) to occupancy bin
%  (this is useful when plotting the tuning variable(s) in units of bins)
%
% CFG:
%
% cfg.binEdges: {1 x nDim} cell array with bin edges (default obtained from
%  min, max in data)
% cfg.occ_dt = 1/30; % time corresponding to each occupancy sample
% cfg.smoothingKernel = [];
% cfg.minOcc = 1; % minimum occupancy (in samples)
%
% MvdM 2014-08-21

% preamble

if ~CheckTSD(tuning_var)
   error('tuning_var input is not a well-formed tsd.'); 
end

cfg_def = [];

nDefaultBins = 100;
cfg_def.nDim = size(tuning_var.data,1); % set up default bins
for iDim = 1:cfg_def.nDim

        mn = min(tuning_var.data(iDim,:));
        mx = max(tuning_var.data(iDim,:));
        
        cfg_def.binEdges{iDim} = linspace(mn,mx,nDefaultBins+1);
end

cfg_def.smoothingKernel = []; % example for 1-D: gausskernel(11,2);
cfg_def.occ_dt = 1/30;
cfg_def.minOcc = 1;
cfg_def.bootstrap = 0;
cfg_def.nBoot = 1000; % number of bootstrap samples to run
cfg_def.bootFrac = 0.9; % fraction of tuning variable data to use for each bootstrap sample
cfg_def.beta0 = 0.5; % initial beta (occupancy) for bins
cfg_def.alpha0 = 2; % expected firing rate (in spks/s)

cfg = ProcessConfig(cfg_def,cfg_in);

% work
switch cfg.nDim
    case 1
        % bin tuning variable
        [tc.occ_hist,tc.pos_idx] = histc(tuning_var.data(1,:),cfg.binEdges{1},2);        
        [tc.occ_hist,tc.pos_idx] = trim_histc(tc.occ_hist,tc.pos_idx);
        
        tc.no_occ_idx = find(tc.occ_hist < cfg.minOcc);
        tc.good_idx = find(tc.occ_hist >= cfg.minOcc);
        
        if ~isempty(cfg.smoothingKernel)
            tc.occ_hist = conv(tc.occ_hist,cfg.smoothingKernel,'same');
        end
        
        tc.occ_hist(tc.no_occ_idx) = NaN;
        tc.occ_hist = tc.occ_hist .* cfg.occ_dt;
        
        % bin spikes
        nCells = length(S.t);
        for iC = 1:nCells
            spk_z = interp1(tuning_var.tvec,tuning_var.data(1,:),S.t{iC},'linear');
            
            if isempty(spk_z), spk_z = zeros(0,1); end % ensures histc generates histogram with all zeros -- MvdM
            
            [spk_hist,~] = histc(spk_z,cfg.binEdges{1},1);
            spk_hist = trim_histc(spk_hist);
            
            if ~isempty(cfg.smoothingKernel)
                spk_hist = conv(spk_hist,cfg.smoothingKernel,'same');
            end
            
            spk_hist(tc.no_occ_idx) = NaN;
            
            tc.tc(iC,:) = spk_hist'./tc.occ_hist;
            
            % if requested, bootstrap
            if cfg.bootstrap
                tdata = tuning_var.data(1,:);
                tt = tuning_var.tvec;
                s = S.t{iC};
                
                s_idx = nearest_idx3(s,tt);
                nT = length(tt);
                
                % intialize bootstrap var
               for iBoot = cfg.nBoot:-1:1
                
                   this_idx = randperm(nT);
                   this_idx = this_idx(1:round(cfg.bootFrac.*nT));
                   this_tdata = tdata(this_idx);
                   this_tt = tt(this_idx);
                   [~,s_toss] = setxor(s_idx,this_idx); % can't use intersect because removes doubles
                   this_s = s; this_s(s_toss) = [];
                   
                   % compute tuning curve over this_tdata and this_s
                   % NOTE this is repeated code, should place in function
                   
                   % bin tuning variable
                   this_occ_hist = histc(this_tdata,cfg.binEdges{1},2);
                   this_occ_hist = trim_histc(this_occ_hist);
                   
                   this_no_occ_idx = find(this_occ_hist < cfg.minOcc);

                   if ~isempty(cfg.smoothingKernel)
                       this_occ_hist = conv(this_occ_hist,cfg.smoothingKernel,'same');
                   end
                   
                   this_occ_hist(this_no_occ_idx) = NaN;
                   this_occ_hist = this_occ_hist .* cfg.occ_dt;
                   
                   % bin spikes
                   this_spk_z = interp1(this_tt,this_tdata,this_s,'linear');
                   
                   if isempty(this_spk_z), this_spk_z = zeros(0,1); end % ensures histc generates histogram with all zeros -- MvdM                
                   
                   [this_spk_hist,~] = histc(this_spk_z,cfg.binEdges{1},1);
                   this_spk_hist = trim_histc(this_spk_hist);
                   
                   if ~isempty(cfg.smoothingKernel)
                       this_spk_hist = conv(this_spk_hist,cfg.smoothingKernel,'same');
                   end
                   
                   this_spk_hist(this_no_occ_idx) = NaN;
                   
                   this_boot_tc(iBoot,:) = this_spk_hist'./this_occ_hist;
                   
               end
               
               % compute SD over bootstrapped TCs
               tc.tcboot(iC,:) = nanstd(this_boot_tc);
            end
        
            
        end
        
        tc.binEdges = cfg.binEdges{1};
        tc.binCenters = cfg.binEdges{1}(1:end-1)+median(diff(cfg.binEdges{1}))/2;
        tc.nBins = length(tc.binEdges)-1;
        
    case 2
        
        % bin tuning variable
        tuning_mat(:,1) = tuning_var.data(1,:)'; % construct input to 2-d histogram
        tuning_mat(:,2) = tuning_var.data(2,:)';

        [tc.beta,~,~,tc.pos_idx] = histcn(tuning_mat,cfg.binEdges{1},cfg.binEdges{2});     
        
        tc.no_occ_idx = find(tc.beta < cfg.minOcc);
        tc.good_idx = find(tc.beta >= cfg.minOcc);

        tc.beta = tc.beta + cfg.beta0;
        
        if ~isempty(cfg.smoothingKernel)
            tc.beta = conv2(tc.beta,cfg.smoothingKernel,'same');
        end

        tc.beta(tc.no_occ_idx) = NaN;
        tc.beta = tc.beta .* cfg.occ_dt;

        % bin spikes
        nCells = length(S.t);
        for iC = 1:nCells

            spk_1 = interp1(tuning_var.tvec,tuning_var.data(1,:),S.t{iC},'linear');
            spk_2 = interp1(tuning_var.tvec,tuning_var.data(2,:),S.t{iC},'linear');
            
            clear spk_mat;
            spk_mat(:,1) = spk_1; spk_mat(:,2) = spk_2;
            spk_hist = histcn(spk_mat,cfg.binEdges{1},cfg.binEdges{2});
            
            if ~isempty(cfg.smoothingKernel)
                spk_hist = conv2(spk_hist,cfg.smoothingKernel,'same');
            end

            spk_hist(tc.no_occ_idx) = NaN;
            
            min_spk = cfg.alpha0 .* cfg.beta0; % expected number of spikes given minimum occupancy
            tc.alpha(iC,:,:) = spk_hist + min_spk;

        end

        tc.binEdges = cfg.binEdges;
        for iDim = 1:cfg.nDim
            tc.binCenters{iDim} = cfg.binEdges{iDim}(1:end-1)+median(diff(cfg.binEdges{iDim}))/2;
            tc.nBins{iDim} = length(tc.binEdges{iDim})-1;
        end

    otherwise
        error('More than 2 tuning dimensions is not yet implemented.');
end




