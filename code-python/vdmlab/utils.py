import numpy as np


def find_nearest_idx(array, val):
    """Finds nearest index in array to value.

    Parameters
    ----------
    array : numpy array
    val : float

    Returns
    -------
    Index into array that is closest to val

    """
    try:
       return (np.abs(array-val)).argmin()
    except:
        import pdb; pdb.set_trace()


def time_slice(spikes, t_start, t_stop):
    """Slices into spike list between a start and stop time.

    Parameters
    ----------
    spikes : list of lists
        Where each inner list contains the spike times for an
        individual neuron. And len(spikes) is the total number
        of neurons.
    t_start : int
        If None, takes the slice from the beginning of the spike
        times.
    t_stop : int
        If None, takes the slice from the end of the spike times.

    Returns
    -------
    sliced_spikes : list of lists
        Where each inner list contains the spike times of interest
        for an individual neuron.

    Raises
    ------
    AssertionError
    When that len(spikes) != len(sliced_spikes) (eg. the number
    of neurons stays the same.
    """
    if t_start is None:
        t_start = -np.inf
    if t_stop is None:
        t_stop = np.inf

    sliced_spikes = []
    for neuron_spikes in spikes:
        indices = (neuron_spikes >= t_start) & (neuron_spikes <= t_stop)
        sliced_spikes.append(neuron_spikes[indices])

    assert(len(spikes) == len(sliced_spikes))

    return sliced_spikes


def idx_in_pos(position, index):
    """ Indexes into position data.

    Keeps x, y, time consistent.

    Parameters
    ----------
    position : dict
        With x, y, time keys
    index : int

    Returns
    -------
    pos : dict
        With x, y, time keys

    """
    pos = dict()
    pos['x'] = position['x'][index]
    pos['y'] = position['y'][index]
    pos['time'] = position['time'][index]
    return pos


def get_sort_idx(tuning_curves):
    """ Finds indices to sort neurons by max firing in tuning curve.

    Parameters
    ----------
    tuning_curves : list of lists
        Where each inner list is the tuning curves for an individual
        neuron.

    Returns
    -------
    sorted_idx : list
        List of integers that correspond to the neuron in sorted order.

    """
    tc_max_loc = []
    for i, neuron_tc in enumerate(tuning_curves):
        tc_max_loc.append((i, np.where(neuron_tc == np.max(neuron_tc))[0][0]))
    sorted_by_tc = sorted(tc_max_loc, key=lambda x: x[1])

    sorted_idx = []
    for idx in sorted_by_tc:
        sorted_idx.append(idx[0])

    return sorted_idx
