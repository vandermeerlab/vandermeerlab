function [ spikes_filtered ] = spikes_filtered_shortcut(spikes, pos_tsd, z, z_iv, expkeys)
% Restricting spikes to Phase 3, along specified trajectory, spiking of
% specified range [min_num_spike <= spikes <= max_num_spike], while rat is
% moving [speed >= threshold]. Results in spikes_running.
%
% Example usage:
% - see shortcut_workflow_tc.m for example
%
% Returns: spikes_filtered

spikes = restrict(spikes, expkeys.phase3(1), expkeys.phase3(2));

% Get neurons that spike between 100 and 4Hz in Phase3
min_num_spike = 100;
max_num_spike = (length(z.tvec)*1/30)*4; %multiply by video sampling rate 1/30
spikes_thres = struct;
neuron = {};
labels = {};
for i = 1:length(spikes.t)
    if cellfun('length',spikes.t(i)) > min_num_spike && cellfun('length',spikes.t(i)) <= max_num_spike
        neuron{end+1} = spikes.t{i};
        labels{end+1} = spikes.label{i};
    end
end
[spikes_thres(:).t] = deal(neuron);
[spikes_thres(:).label] = deal(labels);
[spikes_thres(:).cfg] = deal(spikes.cfg);
[spikes_thres(:).usr] = deal(spikes.usr);

% Restrict neurons to portion that rat is along "ideal" linearized track
spikes_trajectory = restrict(spikes_thres, z_iv);

% Get spikes that happen when rat is running (speed threshold above 10)
speed = tsd;
speed.tvec = pos_tsd.tvec;
speed.data = abs(diff(pos_tsd.data));
speed.label = {'speed'};
cfg_speed = []; 
cfg_speed.method = 'raw'; 
cfg_speed.threshold = 10; 
running_iv = TSDtoIV(cfg_speed, speed);
spikes_running = restrict(spikes_trajectory, running_iv);

% Get spikes that are between 100 and 4Hz in Phase3
spikes_firing.cfg = spikes_running.cfg;
spikes_firing.usr = spikes_running.usr;
neuron_counter = 1;
for i = 1:length(spikes_running.t)
    if length(spikes_running.t{i}) >= min_num_spike && length(spikes_running.t{i}) <= max_num_spike
        spikes_firing.label{neuron_counter} = spikes_running.label{i};
        spikes_firing.t{neuron_counter} = spikes_running.t{i};
        neuron_counter = neuron_counter + 1;
    end
end

spikes_filtered = spikes_running;