import scipy.io as sio


def expand_line(start_pt, stop_pt, line, expand_by=6):
    """ Creates buffer zone around a line.

    Parameters
    ----------
    start_pt : Shapely's Point object
    stop_pt : Shapely's Point object
    line : Shapely's LineString object
    expand_by : int
        This sets by how much you wish to expand the line.
        Defaults to 6.

    Returns
    -------
    zone : Shapely's Polygon object

    """
    line_expanded = line.buffer(expand_by)
    zone = start_pt.union(line_expanded).union(stop_pt)
    return zone


def save_spike_position(spike_position, savepath):
    """ Saves spikes by position as a *.mat for use in matlab

    Parameters
    ----------
    spike_position : list of lists
        Where the inner list contains spike times for an individual
        neuron that are localized to when the animal was in a specified
        position.
    savepath : str
        Entire path for saving this data.
        Eg. 'E:\\code\\vandermeerlab\\code-python\\examples\\emily\\output\\matlab\\r066d1_pedestal_spikes.mat

    """
    for key in spike_position:
        for idx in list(range(len(spike_position[key]))):
            spike_position[key][idx] = spike_position[key][idx].reshape((len(spike_position[key][idx]), 1))
    sio.savemat(savepath, {'spike_pos': spike_position})
