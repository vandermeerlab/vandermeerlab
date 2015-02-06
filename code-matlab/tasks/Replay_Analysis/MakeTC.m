function tc_out = MakeTC(cfg_in,S,pos)
%% MAKETC Make Tuning Curves
%   tc_out = MakeTC(cfg_in,S,pos) returns a struct containing the tuning curves of place
%   cells of spiketrain S. Also contains the spike histogram and occupancy as well as
%   characteristics of the tuning curves (i.e., peaks, order in space)
%
%   INPUT:
%       cfg - input cfg parameters
%       S - input spiketrain
%       Pos - input Position vector binned in cm
%
%   OUTPUT:
%       tc - tuning curve object
%
%   CFG OPTIONS:
%    cfg_def.gkSD = 2;
%    cfg_def.gkWidth = 10;
%    cfg_def.conv2Sec = 1;
%    cfg_def.min_fr = 1; % threshold for detecting place fields in Hz
%    cfg_def.max_meanfr = 5; % mean fr to rule out interneurons
%    cfg_def.maxSpd = 5; % max animal velocity in cm/s
%
% youkitan 2014-12-28

%% Parse cfg parameters
cfg_def.binSize = 1; %position bins in cms
cfg_def.gkSD = 3; % standard deviation for gaussian kernel (in cm)
cfg_def.gkWidth = 20; % width of gaussian kernel in cm (size of place field)
cfg_def.conv2Sec = 1; % convert to seconds/bin
cfg_def.min_fr = 1; % threshold for detecting place fields in Hz
cfg_def.max_meanfr = 5; % mean fr to rule out interneurons

cfg = ProcessConfig2(cfg_def,cfg_in);

%% Detect Place cells and get tuning curves

% Check input data
if any(strcmp(pos.label,'z'))
    pos_mat = getd(pos,'z');
elseif any(strcmp(pos.label,'x'))
    error('Currently does not support 2D place cells')
else
    error('Not a valid position vector')    
end

% Make gaussian kernel and occupancy map 
if isfield(cfg,'rundist')
    nBins = cfg.rundist;
else
    nBins = max(pos_mat);
end

kernel = gausskernel(cfg.gkWidth/cfg.binSize,cfg.gkSD/cfg.binSize);
edges = 0.5:nBins+0.5;

% get occupancy
[occ_hist,pos_idx] = histc(pos_mat,edges);
occ_hist = occ_hist(1:end-1); % ignore last edge bin
occ_hist = conv(occ_hist,kernel,'same');
if cfg.conv2Sec
    sampling_rate = median(diff(pos.tvec)); 
    occ_hist = occ_hist .* sampling_rate;
end
    
% get spiking counts and firing rates
nCells = length(S.t);
for iC = 1:nCells
    if isempty(S.t{iC}) % no spikes
       all_tc(iC,:) = 0;
       continue;
    end
    
    spk_z = interp1(pos.tvec,pos_mat,S.t{iC},'linear');
    spk_hist = histc(spk_z,edges,1);
    spk_hist = spk_hist(1:end-1);
    spk_hist_mat(iC,:) = conv(spk_hist,kernel,'same');

    tc = spk_hist_mat(iC,:)'./occ_hist'; %firing rates = (spikes per bin) / (time spent in bin)
    all_tc(iC,:) = tc;
end

tc = all_tc'; % for compatibility with 2d version
spk_hist = spk_hist_mat;

[template_idx,peak_idx,peak_loc] = DetectPlaceCells1D([],tc);

%store all data in a single struct
tc_out.tc = tc;
tc_out.template_idx = template_idx;
tc_out.peak_idx = peak_idx;
tc_out.peak_loc = peak_loc;
tc_out.occ_hist = occ_hist;
tc_out.spk_hist = spk_hist_mat;

end