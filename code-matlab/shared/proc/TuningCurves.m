function tc_out = TuningCurves(cfg_in,S,tuning_var)
%% TUNING CURVES 
%	tc = TuningCurves(cfg,S,tuning_var) computes tuning curves for cells in S relative to
%   variable(s) tuning_var
%
%	INPUTS:
%       S: ts with spike data
%       tuning_var: tsd with nDim x nSamples tuning variable (e.g. position)
%
%	OUTPUTS:
%       tc_out: tc struct with fields:
%           .tc: nCells x nBins tuning curves
%           .occ_hist: 1 x nBins occupancy (tuning_var sample count)
%           .spk_hist: 1 x nBins spike counts 
%           .tc2D: nCells x nXBins x NYBins tuning curves (for 2D tuning_var only)
%           .binEdges: idx into tc.tc2D specifying which bins are used to make tc.tc
%           .pos_idx = 1 x nSamples idx, assignment of tuning variable(s) to occupancy bin
%          (this is useful when plotting the tuning variable(s) in units of bins)
%
%	CFG:
%       cfg.smoothingKernel = []; % example for 1-D: gausskernel(11,2);
%       cfg.occ_dt = 1/30; % time corresponding to each occupancy sample
%       cfg.minOcc = 1; % minimum occupancy (in samples)
%       cfg.bootstrap = 0; flag, if True use bootstrap to generate tuning curves
%       cfg.nBoot = 1000; % number of bootstrap samples to run
%       cfg.bootFrac = 0.9; % fraction of tuning variable data to use for each bootstrap sample
%
%   see also MakeTC, tc, DetectPlaceCells1D
%   Workflow example: WORKFLOW_PlotOrderedRaster 
%
% MvdM 2014-08-21
% youkitan 2016-12-02 edits: correct help, correct output order, additional features (like MakeTC)
% youkitan edit Feb 2017, restructuring tc format/pipeline

%% cfg preamble and error check

if ~CheckTSD(tuning_var)
   error('tuning_var input is not a well-formed tsd.'); 
elseif ~CheckTS(S)
    error('S is not a well-formed ts.')
elseif size(tuning_var.data,1) > 2
    error('Currently only tuning variables with up to two dimensions are supported.')
end

cfg_def.nBins = 100; %default number of bins unless either specified or using predefined bins

cfg_def.nDim = size(tuning_var.data,1); % set up default bins
for iDim = 1:cfg_def.nDim

        mn = min(tuning_var.data(iDim,:));
        mx = max(tuning_var.data(iDim,:));
        
        cfg_def.binEdges{iDim} = linspace(mn,mx,cfg_def.nBins+1);
end

cfg_def.smoothingKernel = []; % example for 1-D: gausskernel(11,2);
cfg_def.occ_dt = 1/30; % time corresponding to each occupancy sample
cfg_def.minOcc = 1; % minimum occupancy (in samples)
cfg_def.bootstrap = 0;
cfg_def.nBoot = 1000; % number of bootstrap samples to run
cfg_def.bootFrac = 0.9; % fraction of tuning variable data to use for each bootstrap sample

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

%% main body
switch cfg.nDim
    case 1
        % bin tuning variable
        [occ_hist,pos_idx] = histc(tuning_var.data(1,:),cfg.binEdges{1},2);        
        [occ_hist,pos_idx] = trim_histc(occ_hist,pos_idx);
        
        no_occ_idx = find(occ_hist < cfg.minOcc);
        good_idx = find(occ_hist >= cfg.minOcc);
        
        if ~isempty(cfg.smoothingKernel)
            occ_hist = conv(occ_hist,cfg.smoothingKernel,'same');
        end
        
        occ_hist(no_occ_idx) = NaN;
        occ_hist = occ_hist .* cfg.occ_dt;
        
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
            
            spk_hist(no_occ_idx) = NaN;
            spk_hist = spk_hist'; % this ensures that it is the same orientation as occ_hist
            
            tc(iC,:) = spk_hist./occ_hist;
            
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
               tcboot(iC,:) = nanstd(this_boot_tc);
            end
        
            
        end
        
        binEdges = cfg.binEdges{1};
        nBins = length(binEdges)-1;
        binCenters = cfg.binEdges{1}(1:end-1)+median(diff(cfg.binEdges{1}))/2;
        
    case 2
        
        % bin tuning variable
        tuning_mat(:,1) = tuning_var.data(1,:)'; % construct input to 2-d histogram
        tuning_mat(:,2) = tuning_var.data(2,:)';

        [occ_hist,~,~,pos_idx] = histcn(tuning_mat,cfg.binEdges{1},cfg.binEdges{2});     
        
        no_occ_idx = find(occ_hist < cfg.minOcc);
        good_idx = find(occ_hist >= cfg.minOcc);

        if ~isempty(cfg.smoothingKernel)
            occ_hist = conv2(occ_hist,cfg.smoothingKernel,'same');
        end

        occ_hist(no_occ_idx) = NaN;
        occ_hist = occ_hist .* cfg.occ_dt;

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

            spk_hist(no_occ_idx) = NaN;
            
            tc2D(iC,:,:) = spk_hist./occ_hist;
            
            temp = tc2D(iC,:,:);
            tc(iC,:) = temp(good_idx);

        end

        binEdges = cfg.binEdges(1:2);
        for iDim = 1:cfg.nDim
            binCenters{iDim} = cfg.binEdges{iDim}(1:end-1)+median(diff(cfg.binEdges{iDim}))/2;
            nBins{iDim} = length(binEdges{iDim})-1;
        end

    otherwise
        error('More than 2 tuning dimensions is not yet implemented.');
end

%% housekeeping
%tc_out = tc;

% main data
tc_out.tc = tc;
if cfg.bootstrap; tc_out.tcboot = tcboot; end
if exist('tc2D','var'); tc_out.tc2D = tc2D; end
tc_out.occ_hist = occ_hist;
tc_out.spk_hist = spk_hist;

% usr data
tc_out.usr.good_idx = good_idx;
tc_out.usr.no_occ_idx = no_occ_idx;
tc_out.usr.pos_idx = pos_idx;
tc_out.usr.binEdges = binEdges;
tc_out.usr.binCenters = binCenters;
tc_out.usr.nBins = nBins;

% History
tc_out = History(tc_out,mfun,cfg);

end
