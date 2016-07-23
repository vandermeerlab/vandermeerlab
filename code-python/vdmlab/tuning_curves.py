import numpy as np
from scipy import signal
from shapely.geometry import Point

from .utils import find_nearest_idx


def linear_trajectory(pos, ideal_path, trial_start, trial_stop):
    """ Projects 2D positions into an 'ideal' linear trajectory.

    Parameters
    ----------
    pos : dict
        With x, y, time as keys. 2D position information.
    ideal_path : Shapely's LineString object
    trial_start : float
    trial_stop : float

    Returns
    -------
    z : dict
        With position, time as keys

    """
    t_start_idx = find_nearest_idx(np.array(pos['time']), trial_start)
    t_end_idx = find_nearest_idx(np.array(pos['time']), trial_stop)

    pos_trial = dict()
    pos_trial['x'] = pos['x'][t_start_idx:t_end_idx]
    pos_trial['y'] = pos['y'][t_start_idx:t_end_idx]
    pos_trial['time'] = pos['time'][t_start_idx:t_end_idx]

    z = dict(position=[])
    z['time'] = pos_trial['time']
    for point in list(range(len(pos_trial['x']))):
        position = Point(pos_trial['x'][point], pos_trial['y'][point])
        z['position'].append(ideal_path.project(position))
    return z


def tuning_curve(position_z, spike_times, num_bins=100, sampling_rate=1/30.0, filter_type='gaussian', gaussian_std=18):
    """ Computes tuning curves for neurons relative to linear position.

    Parameters
    ----------
    position_z : dict
        With position, time as keys
    spike_times : list of lists
        Where each inner list contains the spike times for an
        individual neuron.
    num_bins : int
        Defaults to 100 if not specified
    sampling_rate : float
        Defaults to 1/30. if not specified.
    filter_type : str, optional
        Defaults to 'gaussian' to smooth with a gaussian filter.
        No smoothing if None.
    gaussian_std : int
        Defaults to 18.

    Returns
    -------
    out_tc : list of lists
        Where each inner list contains the tuning curves for an
        individual neuron.

    Notes
    -----
    Input position_z and spike_times should be from the same time
    period. Eg. when the animal is running on the track.

    """
    linear_start = np.min(position_z['position'])
    linear_stop = np.max(position_z['position'])
    bin_edges = np.linspace(linear_start, linear_stop, num=num_bins)
    bin_centers = np.array((bin_edges[1:] + bin_edges[:-1]) / 2.)
    tc = []
    occupancy = np.zeros(len(bin_centers))
    for position in position_z['position']:
        pos_idx = find_nearest_idx(bin_centers, position)
        occupancy[pos_idx] += sampling_rate
    occupied_idx = occupancy > 0
    for neuron_spikes in spike_times:
        spike_z = np.zeros(len(bin_centers))
        for spike_time in neuron_spikes:
            bin_idx = find_nearest_idx(position_z['time'], spike_time)
            which_bin = find_nearest_idx(bin_centers, position_z['position'][bin_idx])
            spike_z[which_bin] += 1
        firing_rate = np.zeros(len(bin_centers))
        firing_rate[occupied_idx] = spike_z[occupied_idx] / occupancy[occupied_idx]
        tc.append(firing_rate)

    if filter_type == 'gaussian':
        out_tc = []
        filtering_val = 6
        # Normalizing gaussian filter
        gaussian_filter = signal.get_window(('gaussian', gaussian_std), gaussian_std)
        normalized_gaussian = gaussian_filter / np.sum(gaussian_filter)
        for firing_rate in tc:
            out_tc.append(np.convolve(firing_rate, normalized_gaussian, mode='same'))
    else:
        print('Tuning curve with no filter.')
        out_tc = tc

    return out_tc
