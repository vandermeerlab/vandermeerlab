import numpy as np

from .utils import find_nearest_idx


def get_field_idx(neuron, hz_thres, field_dist):
    """ Finds field index from neuron's tuning curve.

    Parameters
    ----------
    neuron : list of floats
        This is the neuron's tuning curve
    hz_thres : int or float
        Everything above this threshold will be included to count as
        part of a given place field.
    field_dist : int or float
        Based on tuning curve bin size. Everything above this distance
        will count as a separate field for a given neuron.

    Returns
    -------
    all_fields : list of list of ints
        Each inner list represents a unique place field, and contains
        the indices where the field is located (in 1D space).

    """
    field_idx = np.where(neuron > hz_thres)[0]

    current_field = []
    all_fields = []
    for idx in range(len(field_idx)):
        if field_idx[idx] < np.max(field_idx):
            if field_idx[idx] not in current_field:
                current_field.append(field_idx[idx])
            if np.abs(field_idx[idx] - field_idx[idx+1]) < field_dist:
                current_field.append(field_idx[idx+1])
            else:
                all_fields.append(current_field)
                current_field = []
    all_fields.append(current_field)
    return all_fields


def find_fields(tuning_curves, hz_thres=5, field_dist=2, max_mean_firing=10):
    """ Finds place field indices for each neuron.

    Parameters
    ----------
    tuning_curves : list of lists of floats
        Where the inner list contains the tuning curves for individual neurons.
    hz_thres : int or float
        Anything above this threshold is considered for being a bin that is part of
        a place field. The default is set at 5.
    field_dist : int or float
        Place fields that are separated by this amount are determined to be separate
        fields. The default is set at 2.
    max_mean_firing : int or float
        Only neurons with a mean firing rate less than this amount are considered for
        having place fields. The default is set to 10.

    Returns
    -------
    all_fields: list of lists of lists of ints
        The outer list contains everything. The first inner list is the neuron, the
        second inner list contains the place field indices. In this way, we can
        represent multiple place fields for individual neurons.
        Eg. 3 neurons that have 2, 1, and 3 place fields respectively would be:
        [[[field], [field]], [[field]], [[field], [field], [field]]]]

    """
    above_thres = []
    for neuron_tc in tuning_curves:
        if np.any(neuron_tc > hz_thres) and (np.mean(neuron_tc < max_mean_firing)):
            above_thres.append(neuron_tc)
        else:
            above_thres.append([])

    all_fields = dict()
    for neuron in range(len(above_thres)):
        if len(above_thres[neuron]) >= 1:
            sorted_fields = get_field_idx(above_thres[neuron], hz_thres, field_dist)
            if len(sorted_fields) >= 1:
                all_fields[neuron] = sorted_fields
    return all_fields


def unique_fields(fields, fields_compare1, fields_compare2):
    """ Finds neurons with and indices of unique fields.

    Parameters
    ----------
    fields : list of lists of lists of ints
        The outer list contains everything. The first inner list is the neuron, the
        second inner list contains the place field indices. In this way, we can
        represent multiple place fields for individual neurons.
        Eg. 3 neurons that have 2, 1, and 3 place fields respectively would be:
        [[[field], [field]], [[field]], [[field], [field], [field]]]]
    fields_compare1 : list of lists of lists of ints
        Same as field, but determined with these neurons on a different tuning curve.
        May also have used more restrictive parameters (eg. hz_thres = 3)
    fields_compare2 : list of lists of lists of ints
        Same as field, but determined with these neurons on a different tuning curve.
        May also have used more restrictive parameters (eg. hz_thres = 3)

    Returns
    -------
    fields_unique : dict of list of lists
        Where the key is the neuron number and the value is a list of place field
        lists.
        Eg. Neurons 7, 3, 11 that have 2, 1, and 3 place fields respectively would be:
        {7: [[field], [field]], 3: [[field]], 11: [[field], [field], [field]]}

    """
    fields_unique = dict()
    for neuron in fields.keys():
        if neuron not in fields_compare1.keys() and neuron not in fields_compare2.keys():
            fields_unique[neuron] = fields[neuron]
    return fields_unique


def sized_fields(fields, min_length=3, max_length=25):
    """ Finds neurons with and indices of properly sized fields.

    Parameters
    ----------
    fields : dict of list of lists
        Where the key is the neuron number and the value is a list of place field
        lists.
        Eg. Neurons 7, 3, 11 that have 2, 1, and 3 place fields respectively would be:
        {7: [[field], [field]], 3: [[field]], 11: [[field], [field], [field]]}
    min_length : int or float
        Based on tuning curve bin size, any fields with fewer than this number of
        bins is not considered a place field.
    max_length : int or float
        Based on tuning curve bin size, any fields with more than this number of
        bins is not considered a place field.

    Returns
    -------
    sized_fields : dict of list of lists
        Where the key is the neuron number and the value is a list of place field
        lists.
        Eg. Neurons 7, 3, 11 that have 2, 1, and 3 place fields respectively would be:
        {7: [[field], [field]], 3: [[field]], 11: [[field], [field], [field]]}

    """
    sized_fields = dict()
    for neuron in fields:
        if any(min_length <= len(field) <= max_length for field in fields[neuron]):
            sized_fields[neuron] = fields[neuron]
    return sized_fields


def get_single_field(fields):
    """ Finds neurons with and indices of single fields.

    Parameters
    ----------
    fields : dict of list of lists
        Where the key is the neuron number and the value is a list of place field
        lists.
        Eg. Neurons 7, 3, 11 that have 2, 1, and 3 place fields respectively would be:
        {7: [[field], [field]], 3: [[field]], 11: [[field], [field], [field]]}

    Returns
    -------
    fields : dict of list of lists
        Where the key is the neuron number and the value is a list of place field
        lists.
        Eg. For the above input, only neuron 3 would be output in this dict:
        {3: [[field]]}

    """
    fields_single = dict()
    for neuron in fields.keys():
        if len(fields[neuron]) == 1:
            fields_single[neuron] = fields[neuron]
    return fields_single


def get_heatmaps(neuron_list, spikes, pos, num_bins=100):
    """ Gets the 2D heatmaps for firing of a given set of neurons.

    Parameters
    ----------
    neuron_list : list of ints
        These will be the indices into the full list of neuron spike times
    spikes : dict
        With times(float), labels (str) as keys
    pos : dict
        With time(float), x(float), y(float) as keys
    num_bins : int
        This will specify how the 2D space is broken up, the greater the number
        the more specific the heatmap will be. The default is set at 100.

    Returns
    -------
    heatmaps : dict of lists
        Where the key is the neuron number and the value is the heatmap for
        that individual neuron.

    """
    xedges = np.linspace(np.min(pos['x'])-2, np.max(pos['x'])+2, num_bins+1)
    yedges = np.linspace(np.min(pos['y'])-2, np.max(pos['y'])+2, num_bins+1)

    heatmaps = dict()
    count = 1
    for neuron in neuron_list:
        field_x = []
        field_y = []
        for spike in spikes['time'][neuron]:
            spike_idx = find_nearest_idx(pos['time'], spike)
            field_x.append(pos['x'][spike_idx])
            field_y.append(pos['y'][spike_idx])
            heatmap, out_xedges, out_yedges = np.histogram2d(field_x, field_y, bins=[xedges, yedges])
        heatmaps[neuron] = heatmap.T
        print(str(neuron) + ' of ' + str(len(neuron_list)))
        count += 1
    return heatmaps

