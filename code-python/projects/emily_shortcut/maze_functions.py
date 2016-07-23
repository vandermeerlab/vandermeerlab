import numpy as np
from shapely.geometry import Point, LineString

import vdmlab as vdm


def get_trial_idx(low_priority, mid_priority, high_priority, feeder1_times, feeder2_times, phase_start, phase_stop):
    start_trials = []
    stop_trials = []

    high_priority_time = []
    mid_priority_time = []
    low_priority_time = []

    f1_idx = 0
    f2_idx = 0

    while f1_idx < len(feeder1_times) and f2_idx < len(feeder2_times):
        if f1_idx == len(feeder1_times):
            start_trial = feeder2_times[f2_idx]
            stop_trial = phase_stop
        elif f2_idx == len(feeder2_times):
            start_trial = feeder1_times[f1_idx]
            stop_trial = phase_stop
        else:
            start_trial = min(feeder1_times[f1_idx], feeder2_times[f2_idx])
            if start_trial in feeder1_times:
                f1_idx += 1
                stop_trial = feeder2_times[f2_idx]
            elif start_trial in feeder2_times:
                f2_idx += 1
                stop_trial = feeder1_times[f1_idx]
        start_trials.append(start_trial)
        stop_trials.append(stop_trial)

        for element in high_priority:
            if np.logical_and(start_trial <= element, element < stop_trial):
                high_priority_time.append(start_trial)
                break
        if start_trial not in high_priority_time:
            for element in mid_priority:
                if np.logical_and(start_trial <= element, element < stop_trial):
                    mid_priority_time.append(start_trial)
                    break
        if start_trial not in high_priority_time and start_trial not in mid_priority_time:
            for element in low_priority:
                if np.logical_and(start_trial <= element, element < stop_trial):
                    low_priority_time.append(start_trial)
                    break

    high_priority_trials = []
    mid_priority_trials = []
    low_priority_trials = []

    for trial in high_priority_time:
        high_priority_trials.append((vdm.find_nearest_idx(np.array(start_trials), trial), 'novel'))
    for trial in mid_priority_time:
        mid_priority_trials.append((vdm.find_nearest_idx(np.array(start_trials), trial), 'shortcut'))
    for trial in low_priority_time:
        low_priority_trials.append((vdm.find_nearest_idx(np.array(start_trials), trial), 'u'))

    trials_idx = dict()
    trials_idx['novel'] = high_priority_trials
    trials_idx['shortcut'] = mid_priority_trials
    trials_idx['u'] = low_priority_trials
    trials_idx['start_trials'] = start_trials
    trials_idx['stop_trials'] = stop_trials

    return trials_idx


def spikes_by_position(spikes, zone, pos_time, pos_x, pos_y):
    spike_position = dict(u=[], shortcut=[], novel=[], other=[])
    counter = 0
    for neuron in spikes:
        neuron_spikes = dict(pedestal=[], u=[], shortcut=[], novel=[], other=[])
        for spike in neuron:
            pos_idx = vdm.find_nearest_idx(pos_time, spike)
            point = Point([pos_x[pos_idx], pos_y[pos_idx]])
            if zone['pedestal'].contains(point) or zone['uped'].contains(point) or zone['shortped'].contains(point) or zone['novelped'].contains(point):
                neuron_spikes['pedestal'].append(np.asscalar(spike))
                continue
            elif zone['u'].contains(point) or zone['ushort'].contains(point) or zone['unovel'].contains(point):
                neuron_spikes['u'].append(np.asscalar(spike))
                continue
            elif zone['shortcut'].contains(point):
                neuron_spikes['shortcut'].append(np.asscalar(spike))
                continue
            elif zone['novel'].contains(point):
                neuron_spikes['novel'].append(np.asscalar(spike))
                continue
            else:
                neuron_spikes['other'].append(np.asscalar(spike))
        spike_position['u'].append(np.array(neuron_spikes['u']))
        spike_position['shortcut'].append(np.array(neuron_spikes['shortcut']))
        spike_position['novel'].append(np.array(neuron_spikes['novel']))
        spike_position['other'].append(np.array(neuron_spikes['other']))
        counter += 1
        print(str(counter) + ' of ' + str(len(spikes)) + ' neurons completed!')
    return spike_position


def get_zones(info, t_start, t_stop):
    pos = info.get_pos(info.pxl_to_cm)
    # Slicing position to only Phase 3
    t_start_idx = vdm.find_nearest_idx(np.array(pos['time']), t_start)
    t_end_idx = vdm.find_nearest_idx(np.array(pos['time']), t_stop)

    sliced_pos = dict()
    sliced_pos['x'] = pos['x'][t_start_idx:t_end_idx]
    sliced_pos['y'] = pos['y'][t_start_idx:t_end_idx]
    sliced_pos['time'] = pos['time'][t_start_idx:t_end_idx]

    # Here I define the ideal trajectories in cm that I project onto
    # to make the 1D linear position.
    u_line = LineString(info.u_trajectory)
    shortcut_line = LineString(info.shortcut_trajectory)
    novel_line = LineString(info.novel_trajectory)

    u_start = Point(info.u_trajectory[0])
    u_stop = Point(info.u_trajectory[-1])
    shortcut_start = Point(info.shortcut_trajectory[0])
    shortcut_stop = Point(info.shortcut_trajectory[-1])
    novel_start = Point(info.novel_trajectory[0])
    novel_stop = Point(info.novel_trajectory[-1])

    zones = dict()
    zones['u'] = vdm.expand_line(u_start, u_stop, u_line)
    zones['shortcut'] = vdm.expand_line(shortcut_start, shortcut_stop, shortcut_line)
    zones['novel'] = vdm.expand_line(novel_start, novel_stop, novel_line)
    zones['ushort'] = zones['u'].intersection(zones['shortcut'])
    zones['unovel'] = zones['u'].intersection(zones['novel'])

    u_idx = []
    shortcut_idx = []
    novel_idx = []
    other_idx = []
    for pos_idx in list(range(len(sliced_pos['time']))):
        point = Point([sliced_pos['x'][pos_idx], sliced_pos['y'][pos_idx]])
        if zones['u'].contains(point) or zones['ushort'].contains(point) or zones['unovel'].contains(point):
            u_idx.append(pos_idx)
        elif zones['shortcut'].contains(point):
            shortcut_idx.append(pos_idx)
        elif zones['novel'].contains(point):
            novel_idx.append(pos_idx)
        else:
            other_idx.append(pos_idx)

    spike_pos = dict()
    spike_pos['u'] = vdm.idx_in_pos(sliced_pos, u_idx)
    spike_pos['shortcut'] = vdm.idx_in_pos(sliced_pos, shortcut_idx)
    spike_pos['novel'] = vdm.idx_in_pos(sliced_pos, novel_idx)
    spike_pos['other'] = vdm.idx_in_pos(sliced_pos, other_idx)

    return spike_pos
