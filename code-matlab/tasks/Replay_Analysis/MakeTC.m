function tc_out = MakeTC(cfg_in,S,pos)
%% MAKETC Make Tuning Curves
%   tc_out = MakeTC(cfg_in,S,pos) returns a struct containing the tuning curves of place
%   cells of spiketrain S. Also contains the spike histogram and occupancy as well as
%   characteristics of the tuning curves (i.e., peaks, order in space)
%
%   INPUT:
%       cfg - input cfg parameters
%       S - input spiketrain
%       pos - input position vector (assumed in cm)
%
%   OUTPUT:
%       tc - tuning curve object
%
%   CFG OPTIONS:
%    cfg_def.gkSD = 2;
%    cfg_def.gkWidth = 10;
%    cfg_def.conv2Sec = 1;
%    cfg_def.p_thr = 5; % threshold for detecting place fields in Hz
%    cfg_def.max_meanfr = 5; % mean fr to rule out interneurons
%    cfg_def.verbose = 1; 1 display command window text, 0 don't
%
% youkitan 2014-12-28, MvdM edits
% youkitan 2016-12-02 youkitan edit: tc size check

%% Parse cfg parameters
cfg_def.binSize = 1; % how many cm in a bin
cfg_def.debug = 0;
cfg_def.gkSD = 3; % standard deviation for gaussian kernel (in cm)
cfg_def.gkWidth = 25; % width of gaussian kernel window (in cm)
cfg_def.conv2Sec = 1; % convert to seconds/bin
cfg_def.p_thr = 5; % min threshold for detecting candidate place fields in Hz
cfg_def.max_meanfr = 10; % mean fr to rule out interneurons
cfg_def.minSize = 4; % minimum size of field (in bins)
%cfg_def.maxSize = 20; % maximum size of field (in bins)
cfg_def.edges = []; % edge vector to use for binning; defaults to min(pos):binSize:max(pos)
cfg_def.nSpikesInField = 20; % minimum number of spikes in field
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose; disp([mfun,': gathering tuning curve and place cell data...']); end

%% Detect place cells and get tuning curves

% Check input data
if any(strcmp(pos.label,'z'))
    pos_mat = getd(pos,'z');
elseif any(strcmp(pos.label,'x'))
    error('Currently does not support 2D place cells')
else
    error('Not a valid position vector')    
end

% define bin edges
if isempty(cfg.edges)
    cfg.edges = min(pos_mat)-0.5:cfg.binSize:max(pos_mat)+0.5;
end
nBins = length(cfg.edges)-1;

% smoothing kernel
kernel = gausskernel(cfg.gkWidth/cfg.binSize,cfg.gkSD/cfg.binSize);

% get occupancy
[occ_hist,pos_idx] = histc(pos_mat,cfg.edges);
occ_hist = trim_histc(occ_hist); % ignore last edge bin
occ_hist = conv(occ_hist,kernel,'same');
if cfg.conv2Sec
    Fs = median(diff(pos.tvec));
    fprintf('MakeTC.m: video Fs %.2f detected.\n',1./Fs);
    occ_hist = occ_hist .* Fs;
end
    
% get spiking counts and firing rates
nCells = length(S.t);
all_tc = nan(nCells,nBins);
spk_hist_mat = nan(nCells,nBins);

for iC = 1:nCells
    if isempty(S.t{iC}) % no spikes
       all_tc(iC,:) = 0;
       spk_hist_mat(iC,:) = 0; %this ensures that both matrices are the same size!
       continue;
    end
    
    spk_z = interp1(pos.tvec,pos_mat,S.t{iC},'linear');
    spk_hist = histc(spk_z,cfg.edges,1);
    spk_hist = spk_hist(1:end-1);
    spk_hist_mat(iC,:) = conv(spk_hist,kernel,'same');

    tc = spk_hist_mat(iC,:)'./occ_hist'; %firing rates = (spikes per bin) / (time spent in bin)
    all_tc(iC,:) = tc;
end

tc_temp.tc = all_tc;
tc_temp.usr = S.usr;
spk_hist = spk_hist_mat;

% detect place fields
detect_cfg = [];
detect_cfg.S = S;
detect_cfg.pos = pos;
detect_cfg.nSpikesInField = cfg.nSpikesInField;
%detect_cfg.p_thr = cfg.p_thr; % threshold for detecting place fields in Hz
detect_cfg.max_meanfr = cfg.max_meanfr;
detect_cfg.minSize = cfg.minSize;
%detect_cfg.maxSize = cfg.maxSize; 
detect_cfg.debug = cfg.debug;
tc_out = DetectPlaceCells1D(detect_cfg,tc_temp);

template_idx = tc_out.template_idx;
peak_idx = tc_out.peak_idx;
peak_loc = tc_out.peak_loc;

% create field template
fpeak_val = []; fpeak_idx = [];
iPeak = 1;
for iP = 1:length(peak_loc)
    fpeak_val(iPeak:iPeak+length(peak_loc{iP})-1) = peak_loc{iP};
    fpeak_idx(iPeak:iPeak+length(peak_loc{iP})-1) = template_idx(iP);
    
    iPeak = iPeak + length(peak_loc{iP});
end
[fpeak_val,sort_idx] = sort(fpeak_val,'ascend');
fpeak_idx = fpeak_idx(sort_idx);

%error check before finishing
if any(size(tc_out.tc') ~= size(spk_hist_mat))
    error('the output does not match!')
end


%store all data in a single struct
tc_out.tc = tc_out.tc';
%tc_out.template_idx = template_idx;
tc_out.field_template_idx = fpeak_idx;
tc_out.field_loc = fpeak_val;
tc_out.peak_idx = peak_idx;
tc_out.peak_loc = peak_loc;
tc_out.occ_hist = occ_hist;
tc_out.spk_hist = spk_hist_mat;

end