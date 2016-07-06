function data_out = AMPX_get_phase(cfg_in, data_ft, ExpKeys)
% function AMPX_get_phase (cfg_in, data_ft, ExpKeys)
% Takes in trialified data and extracts the phase using the Matlab cpsd
% method.  cpsd is used instead of the Fieldtrip funciton because it is
% more precise.  This will return a structure with the phase difference
% for a given set of references across all the trials.
%
%To Do Prior to this function:
%   ensure that the ExpKeys is complete and contains the channel layout and
%   the bad channels.  The Expkeys script will be called and used
%   throughout this function.
%
%INPUTS:
%   trial_data: trialified data in the fieldtrip format (trial_data).  This will include
%   a header, tvec, trial data (nchan x length)
%
%
% OUTPUTS:
%   phase_distrib = [1 x nTrials x reference channels] cell array with each element being the
%   phase for each channel
%   phase_distrib_avg = [8 x 8] matrix of the average
%% parameters
cfg_def = [];
cfg_def.spec_window = 128;
cfg_def.NFFT = 1024;
cfg_def.ref = Naris_BestChan_remap(ExpKeys, 'location', 'vl');
cfg_def.ref_loc = 'vl'; % which corner of the probe to use as a reference (can be 'vl' (ventrolateral), 'vm' (ventromedial), 'dl' (dorsolateral), 'dm' (dorsal medial)
cfg_def.output = 'rad'; % output in radians ('rad') or degrees ('deg')
cfg = ProcessConfig(cfg_def, cfg_in);
%% determine the reference electrode

% [ref,ref_loc] = Naris_BestChan_remap(ExpKeys);
% ref = ref_loc.(cfg.ref_loc);

%% create some empty cells
% all_coh_spec = cell(1,length(data_ft.trial)); % will hold the phase for each channel relative to the reference.


%% compute the phase difference between each channel and the reference using cross-spectral power desnity.
% a negative values represent lags, while positives are leads.
for itrial = 1:length(data_ft.trial)
    for ii = 64:-1:1
        [Cxy,F] = cpsd(data_ft.trial{itrial}(cfg.ref,:),data_ft.trial{itrial}(ii,:),hamming(cfg.spec_window/2),cfg.spec_window/4,cfg.NFFT,data_ft.hdr.Fs);
        coh_spec_fake = -angle(Cxy); %higher value means leading. outputs radians
        all_coh_spec{itrial}(1,ii) = circ_mean(coh_spec_fake(nearest(F,cfg.freq(1)):nearest(F,cfg.freq(2))));
    end
    for bad = ExpKeys.BadChannels;
        for ii = 64:-1:1
            if ismember(ExpKeys.Probe_layout(ii), bad)
                all_coh_spec{itrial}(1,ii) = nan;
            end
        end
    end
    data_out.phase{itrial} = reshape(all_coh_spec{itrial}, 8, 8);
end


%% compute the average
% convert to 1 x 64 for each trial since circ_median can only hand dim of 1
% or 2.
%median method
for itrial = length(data_ft.trial):-1:1
    data_long(itrial,:) = reshape(data_out.phase{itrial},1,64);
end

data_out.phase_avg_median = reshape(circ_median(data_long, 1),8,8);

%circ_mean method
all_trials_3d = cat(3,data_out.phase{:});
data_out.phase_avg_mean = circ_mean(all_trials_3d, [], 3);


