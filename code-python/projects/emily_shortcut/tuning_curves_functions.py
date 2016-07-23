import os
import numpy as np
import pickle
from shapely.geometry import Point, LineString

import vdmlab as vdm

from maze_functions import spikes_by_position


def linearize(info, pos, t_start, t_stop, expand_by=6):
    """Finds linear and zones for ideal trajectories.

        Parameters
        ----------
        info : module
            Contains session-specific information.
        pos : dict
            With x, y, time as keys. Each value is a np.array.
        t_start : float
        t_stop : float
        expand_by : int or float
            This is how much you wish to expand the line to fit
            the animal's actual movements. Default is set to 6.

        Returns
        -------
        linear : dict
            With u, shortcut, novel keys. Each value is a unique
            Shapely LineString object.
        zone : dict
            With 'ushort', 'u', 'novel', 'uped', 'unovel', 'pedestal',
            'novelped', 'shortcut', 'shortped' keys.
            Each value is a unique Shapely Polygon object.

        """
    # Slicing position to only Phase 3
    t_start_idx = vdm.find_nearest_idx(pos['time'], t_start)
    t_end_idx = vdm.find_nearest_idx(pos['time'], t_stop)
    sliced_pos = dict()
    sliced_pos['x'] = pos['x'][t_start_idx:t_end_idx]
    sliced_pos['y'] = pos['y'][t_start_idx:t_end_idx]
    sliced_pos['time'] = pos['time'][t_start_idx:t_end_idx]

    u_line = LineString(info.u_trajectory)
    shortcut_line = LineString(info.shortcut_trajectory)
    novel_line = LineString(info.novel_trajectory)

    u_start = Point(info.u_trajectory[0])
    u_stop = Point(info.u_trajectory[-1])
    shortcut_start = Point(info.shortcut_trajectory[0])
    shortcut_stop = Point(info.shortcut_trajectory[-1])
    novel_start = Point(info.novel_trajectory[0])
    novel_stop = Point(info.novel_trajectory[-1])
    pedestal_center = Point(info.path_pts['pedestal'][0], info.path_pts['pedestal'][1])
    pedestal = pedestal_center.buffer(expand_by*2.2)

    def expand_line(start_pt, stop_pt, line, expand_by):
        line_expanded = line.buffer(expand_by)
        zone = start_pt.union(line_expanded).union(stop_pt)
        return zone

    zone = dict()
    zone['u'] = expand_line(u_start, u_stop, u_line, expand_by)
    zone['shortcut'] = expand_line(shortcut_start, shortcut_stop, shortcut_line, expand_by)
    zone['novel'] = expand_line(novel_start, novel_stop, novel_line, expand_by)
    zone['ushort'] = zone['u'].intersection(zone['shortcut'])
    zone['unovel'] = zone['u'].intersection(zone['novel'])
    zone['uped'] = zone['u'].intersection(pedestal)
    zone['shortped'] = zone['shortcut'].intersection(pedestal)
    zone['novelped'] = zone['novel'].intersection(pedestal)
    zone['pedestal'] = pedestal

    u_idx = []
    shortcut_idx = []
    novel_idx = []
    other_idx = []
    for pos_idx in range(len(sliced_pos['time'])):
        point = Point([sliced_pos['x'][pos_idx], sliced_pos['y'][pos_idx]])
        if zone['u'].contains(point) or zone['ushort'].contains(point) or zone['unovel'].contains(point):
            u_idx.append(pos_idx)
        elif zone['shortcut'].contains(point) or zone['shortped'].contains(point):
            shortcut_idx.append(pos_idx)
        elif zone['novel'].contains(point) or zone['novelped'].contains(point):
            novel_idx.append(pos_idx)
        else:
            other_idx.append(pos_idx)

    u_pos = vdm.idx_in_pos(sliced_pos, u_idx)
    shortcut_pos = vdm.idx_in_pos(sliced_pos, shortcut_idx)
    novel_pos = vdm.idx_in_pos(sliced_pos, novel_idx)
    other_pos = vdm.idx_in_pos(sliced_pos, other_idx)

    linear = dict()
    linear['u'] = vdm.linear_trajectory(u_pos, u_line, t_start, t_stop)
    linear['shortcut'] = vdm.linear_trajectory(shortcut_pos, shortcut_line, t_start, t_stop)
    linear['novel'] = vdm.linear_trajectory(novel_pos, novel_line, t_start, t_stop)

    return linear, zone


def get_tc(info, pos, pickle_filepath):
    """Loads saved tuning curve if it exists, otherwise computes tuning curve.

        Parameters
        ----------
        info : module
            Contains session-specific information.
        pos : dict
            With x, y, time as keys. Each value is a np.array.
        pickle_filepath: str
            Absolute (or relative) location of where tuning_curve.pkl files are saved.

        Returns
        -------
        tc : dict
            With u, shortcut, novel keys. Each value is a list of list, where
            each inner list represents an individual neuron's tuning curve.

        """
    pickled_tc = pickle_filepath + info.session_id + '_tuning_curves_phase3.pkl'
    if os.path.isfile(pickled_tc):
        with open(pickled_tc, 'rb') as fileobj:
            tc = pickle.load(fileobj)
    else:
        t_start = info.task_times['phase3'][0]
        t_stop = info.task_times['phase3'][1]

        spikes = info.get_spikes()

        linear, zone = linearize(info, pos, t_start, t_stop)

        pickled_spike_pos = pickle_filepath + info.session_id + '_spike_position_phase3.pkl'
        if os.path.isfile(pickled_spike_pos):
            with open(pickled_spike_pos, 'rb') as fileobj:
                spike_position = pickle.load(fileobj)
        else:
            sliced_spikes = vdm.time_slice(spikes['time'], t_start, t_stop)
            spike_position = spikes_by_position(sliced_spikes, zone, pos['time'], pos['x'], pos['y'])
            with open(pickled_spike_pos, 'wb') as fileobj:
                pickle.dump(spike_position, fileobj)

        tc = dict()
        tc['u'] = vdm.tuning_curve(linear['u'], spike_position['u'], num_bins=47)
        tc['shortcut'] = vdm.tuning_curve(linear['shortcut'], spike_position['shortcut'], num_bins=47)
        tc['novel'] = vdm.tuning_curve(linear['novel'], spike_position['novel'], num_bins=47)
        with open(pickled_tc, 'wb') as fileobj:
            pickle.dump(tc, fileobj)

    return tc


def get_odd_firing_idx(tuning_curve, max_mean_firing=10):
    """Find indices where neuron is firing too much to be condidered a place cell

    Parameters
    ----------
    tuning_curve :
    max_mean_firing : int or float
        A neuron with a max mean firing above this level is considered to have odd
        firing and it's index will be added to the odd_firing_idx.

    Returns
    -------
    odd_firing_idx : list of ints
        Where each int is an index into the full list of neurons.

        """
    odd_firing_idx = []
    for idx in range(len(tuning_curve)):
        if (np.mean(tuning_curve[idx]) > max_mean_firing):
            odd_firing_idx.append(idx)
    return odd_firing_idx
