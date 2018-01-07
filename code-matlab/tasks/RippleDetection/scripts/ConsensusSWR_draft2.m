%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%             Implement Consensus SWR Detection Method                %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% See full description of method in:
%           Yu et al, 2017 - Distinct hippocampal-cortical memory 
%           representations for experiences associated with movement versus
%           immobility https://elifesciences.org/articles/27621
%
% Shortcut must include the Ripple branch of the vandermeerlab codebase
% (will not work on any other branches)
%
% aacarey Jan 2017 

%% Load LFPs

% Set current directory
cd('D:\data\R064\R064-2015-04-23')

% Load all available LFPs into a [1xn struct]
% Loading all available because "A consensus SWR envelope was calculated by
% taking the median of the envelopes across all available tetrodes. "
% "Only days where at least three tetrodes were in or near the CA1 cell
% layer were used for the analysis."
% R064 has more than 3 tetrodes in the layer
flist = FindFiles('*CSC*');

for iFile = 1:length(flist)
    cfg.fc = {flist{iFile}};
    LFP(iFile) = LoadCSC(cfg);
end

%% Detect SWRs (in TSD form, thresholding to be done later)

% SET CONFIG OPTIONS

% OldWizard (SWR detector relying on Hilbert transform)
% "The raw CA1 LFP was referenced to an electrode in corpus callosum and
% then filtered (150–250 Hz) to isolate the SWR band.  The SWR envelope was
% then obtained using the Hilbert transform and then smoothed by convolving 
% with a Gaussian kernel (? = 4 ms)"
% R064's LFPs are already referenced to CC, so no need to do it here.
cfg_wizard = [];
cfg_wizard.rippleband = [150 250]; % Set ripple band here
cfg_wizard.smooth = 1;
cfg_wizard.kernel = gausskernel(60,8); % at 2 kHz sample rate, 8 samples is 4 ms sigma; the 60 here doesn't really matter, but that's the width of the kernel in samples
cfg_wizard.weightby = 'amplitude'; % 'power' or 'amplitude'; note 'amplitude' is the signal envelope

% Perform the detection and gather the outputs into [1xn structs]
for iLFP = 1:length(flist)
     
    % Hilbs
    %--detect ripples with hilbert transform
    SWRTSD_hilbert(iLFP) = OldWizard(cfg_wizard,LFP(iLFP));
end

%% Get consensus TSD for the detectors

% Combine the data values to get a consensus on the location of SWRs (paper uses median)
% "A consensus SWR envelope was calculated by taking the median of the
% envelopes across all available tetrodes. "
cfg = [];
cfg.method = 'median';
SWRTSD_hilbert_consensus = MergeTSD(cfg,SWRTSD_hilbert);

%% Load some data for plotting
cfg = [];  cfg.load_questionable_cells = 0;
S = LoadSpikes(cfg);
 
LoadExpKeys; cfg = []; cfg.fc = ExpKeys.goodSWR(1);
CSC = LoadCSC(cfg);

% Load human-rated intervals
IV_annotated = loadpop('R064-2015-04-23-manualIVann.mat');
 
% Rescale the TSDs so they will be visible in the plot (too small
% otherwise)
SWR_hilbert = NormalizeTSD([],SWRTSD_hilbert_consensus,'range',[0 length(S.t)]);

%% Plot the detectors along with human-identified intervals

cfg_mr = []; cfg_mr.lfp = CSC; cfg_mr.evt = IV_annotated;
MultiRaster(cfg_mr,S); overplot([],SWR_hilbert);

%% Get ROC data
% ROC can be used to recommend the ideal threshold, as long as there are
% annotated intervals to compare against. This is not feasible for every
% recorded session that undergoes analysis, but in our annotated datasets,
% we can use the ROC functions as a form of validation for the methods
% outline in Yu et al. When a method is validated, it could then be assumed
% to perform at least as well in new data sets. ("safe")

ROCdata_hilbert_consensus = GetROCdata([],IV_annotated,SWRTSD_hilbert_consensus);

disp(' ')
disp('CONSENSUS SWR DETECTION (HILBERT TRANSFORM)')
Score_hilbert_consenses = ScoreROCdata([],ROCdata_hilbert_consensus);

%% Recommend threshold based on noise distribution

% The ROC data functions recommend a threshold to use on the whole session
% to optimize the false positive / false neagtive rate, let's see what
% recommendation the noise distribtution method suggests to use, given the
% whole TSD:
disp(' '); disp('THRESHOLD (WHOLE SESSION)')
threshold_whole_session = NoiseDistributionThreshold([],SWRTSD_hilbert_consensus);

%% Load Position data
% Yu et al restrict the TSD to times when the animal is more quiescent
% before obtaining the noise distribution threshold
% "normalized the consensus envelope power to times during immobility
% (speed <4 cm/s)."

LoadExpKeys
cfg = [];
cfg.convFact = ExpKeys.convFact; % get conversion factor for pixels to centimeters
pos = LoadPos(cfg);

% calculate linear speed
linear_speed = getLinSpd([],pos);

% find intervals where the linear speed is below threshold (4 cm/s in Yu et
% al)
cfg = []; cfg.method = 'raw'; cfg.threshold = 4; % in cm/sec, because of cfg.conFact during pos loading
cfg.operation = '<';
iv_slow = TSDtoIV(cfg,linear_speed); % intervals with speed above thresh

SWR_slow = restrict(SWRTSD_hilbert_consensus,iv_slow);


%% Recommend threshold for data that has been restricted to periods of 
% quiescence
% This is how Yu et al handed in their data to get the threshold
% recommendation

ROCdata_hilbert_consensus_slow = GetROCdata([],IV_annotated,SWR_slow);

disp(' ')
disp('CONSENSUS SWR DETECTION (HILBERT TRANSFORM)')
Score_hilbert_consenses = ScoreROCdata([],ROCdata_hilbert_consensus_slow);

disp(' '); disp('THRESHOLD (PERIODS OF QUIESCENCE)')
threshold_slow = NoiseDistributionThreshold([],SWR_slow);

%% Evaluate the proportion of human-rated events that have been detected, 
% as well as the proportion of total detected events that are false positives

% It seems that they perform threshold recommendation on the consensus signal
% envelope, but perform thresholding on the consensus power envelope:
% "The SWR envelope was then obtained using the Hilbert transform and then 
% smoothed by convolving..."
% Then description of noise distribution calculations
% " SWRs were then detected when the power of the consensus envelope 
% exceeded the detection threshold for a minimum of 20 ms."

% Get power envelope
SWRTSD_hilbert_consensus_power = SWRTSD_hilbert_consensus;
SWRTSD_hilbert_consensus_power.data = SWRTSD_hilbert_consensus_power.data.^2; 

% Now perform thresholding, producing intervals that define the start and
% end times of SWRs
cfg = [];
cfg.method = 'zscore'; cfg.threshold = threshold_slow;
SWRIV_consensus_slow = TSDtoIV2(cfg,SWRTSD_hilbert_consensus_power);

% Exclude intervals that are too short in duration
cfg = []; cfg.mindur = 0.02;
SWRIV_consensus_slow = RemoveIV(cfg,SWRIV_consensus_slow);

% Exclude intervals when the animal is moving quickly
cfg = []; cfg.straddle = 0;
SWRIV_detected = RestrictIV(cfg,SWRIV_consensus_slow,iv_slow);

% Display hit rate results
% 1's are intervals that aacarey were very confident were real SWRs
% 5's are intervals that aacarey thinks could be false positives, however their
% spectrograms show slightly higher power in the ripple band than
% surrounding regions
% The false positive rate is shown as a proportion of the whole set of
% detected intervals. (Are they 5% false positives, 50% false positives,
% etc)
hitrates = HitRate([],IV_annotated,SWRIV_detected); 
