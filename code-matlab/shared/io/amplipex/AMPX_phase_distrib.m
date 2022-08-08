\function [phase] = AMPX_phase_distrib (trial_data, varargin)
% function AMPX_phase_distrib (trial_data, varargin)
% Takes in trialified data and extracts the phase using the Matlab cpsd
% method.  cpsd is used instead of the Fieldtrip funciton because it is
% more reliable.  This will return a structure with the phase difference
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
%   phase_distrib_avg = [8 x 8 x reference channels] matrix of the average
%   phase value for each reference site.
%


%% Define the variables
run(FindFile('*Keys.m'))

min_freq = 45; % low gamma [45 65], high gamma [70 85], high_freq [450 490]
max_freq = 65;
spec_window = 128;
trials_def = 1:length(trial_data.trial);
y_label = -pi:1/10*pi:pi;
ref_chan = ExpKeys.Ref_chan;

extract_varargin;

%% use the cpsd method to get the phase  of the

all_coh_spec = cell(64,1); all_phase_spec_r = cell(64,1); temp_1 = zeros(1,64);
phase.phase_distrib = cell(max(trials_def),1);
phase_distrib_avg = cell(length(ref_chan),1);
[~, all_best] = AMPX_BestChan(ExpKeys);
ref_chan = all_best.vl;
ref = find(ExpKeys.Probe_layout == ref_chan);
disp(['Reference Channel:  ' num2str(ref_chan)])
for itrial = trials_def;
    all_coh_spec{ref} = zeros(64,1);
    for ii = 64:-1:1
        [Cxy,F] = cpsd(trial_data.trial{itrial}(ref_chan,:),trial_data.trial{itrial}(ii,:), rectwin(length(trial_data.trial{itrial}(ii,:))), 0,256*10, trial_data.hdr.Fs); %,hamming(spec_window/2),spec_window/4,1024,trial_data.hdr.Fs);
        coh_spec = -angle(Cxy); %higher value means leading. outputs radians
        all_coh_spec{ref}(ii,1) = circ_mean(coh_spec(nearest(F,min_freq):nearest(F,max_freq)));
        
    end
    temp_1(1,33:64) = all_coh_spec{ref}(33:64,1)'; temp_1(1,1:32) = all_coh_spec{ref}(1:32,1)';
    for bad = ExpKeys.BadChannels;
        for ii = 64:-1:1
            if ExpKeys.Probe_layout(ii) == bad
                temp_1(1,ii) = nan;
            end
        end
    end
    all_phase_spec_r{ref} = reshape(temp_1,8,8);
    phase.phase_distrib{itrial} = all_phase_spec_r;
%     all_phase_spec_r{ref} = reshape(temp_ref, 8,8);
end



ref_phases = [];
for itrial = trials_def
    ref_phases = cat(3, ref_phases, phase.phase_distrib{itrial}{ref});
end

    phase.phase_distrib_avg = circ_mean(ref_phases, [],3);

% for ii = 64:-1:1
%     loop_ctn = 1;
%     for itrial = trials_def
%         ref_phases(loop_ctn) = phase.phase_distrib{itrial}{ref}(ii);
%         loop_ctn = loop_ctn +1;
%     end
%     phase_distrib_avg{ref}(ii) = circ_mean(ref_phases');
% end

%% add the paramters to the power structure
phase.params.freq = [min_freq max_freq];
phase.params.window_size = ((trial_data.cfg.trl(1,2)/trial_data.fsample)-(trial_data.cfg.trl(1,1)/trial_data.fsample))/2;
phase.params.ref_channels = ref_chan;





