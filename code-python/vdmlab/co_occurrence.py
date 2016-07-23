import numpy as np


def spike_counts(spike_times, interval_times, window=None):
    """ Get spike counts for specific interval.

    Parameters
    ----------
    spike_times :
    interval_times : dict
        With start(int or float), stop(int or float) as keys
    window : float
        When window is set, takes the spike times for this window length
        around the center of the interval times. The default is set to None.
        Suggested window size is 0.1.

    Returns
    -------
    count_matrix : np.array
        num_neurons x num_bins

    """
    intervals = np.vstack((interval_times['start'], interval_times['stop']))
    bin_centers = np.mean(intervals, axis=0)

    if window is not None:
        intervals = np.vstack([bin_centers-(window*0.5), bin_centers+(window*0.5)])

    num_neurons = len(spike_times)
    count_matrix = np.zeros((num_neurons, intervals.shape[1]))

    for i, (start, stop) in enumerate(zip(intervals[0], intervals[1])):
        for neuron in range(num_neurons):
            count_matrix[neuron][i] = ((start <= spike_times[neuron]) & (spike_times[neuron] <= stop)).sum()

    return count_matrix


def cooccur_probabilities(count_matrix, shuffle=10000):
    """ Obtain expected and observed co-occurrence probabilities

    Parameters
    ----------
    count_matrix : np.array
        num_neurons x num_bins
    shuffle : int
        Number of shuffles to compare your data against. The default
        is set to 10000.

    Returns
    -------
    output : dict
        Where the keys are p0, p1, p2, p3, p4, p5. P1 through p4 are
        np.arrays of length num_neurons

    Notes
    -----
    p0 : probability (fraction of time bins) each neuron is active.
    p1 : expected co-occurrence under independence assumption
        ..math:: p(x,y) = p(x) * p(y)
    p2 : Observed conditional probability
        ..math:: p(x|y)
    p3 : Observed co-occurrence (joint) probability
        ..math:: p(x,y)
    p4 : z-score of p3 against shuffled data
    p5 : Observed co-occurrence (joint) probability of shuffled data
        ..math:: p(x,y)

        """
    # Boolean for if at least on spike happened in this time bin
    count_matrix[count_matrix > 1] = 1

    # Get expected co-occurrence under independence assumption
    # fraction of bins each cell participates in individually
    p0 = np.mean(count_matrix, axis=1)

    num_neurons = count_matrix.shape[0]
    # Expected co-occurrence, multiply single cell probabilities
    p1 = np.zeros((num_neurons, num_neurons))
    for i in range(num_neurons):
        for j in range(num_neurons):
            p1[i][j] = p0[i] * p0[j]
    # remove probability of cell co-occurring with itself
    p1[np.eye(len(p1), dtype=bool)] = np.nan

    # Observed conditional probabilities
    # (probability of cell y active given that cell x is active)
    p2 = np.zeros((num_neurons, num_neurons))
    for i in range(num_neurons):
        col_idx = count_matrix[i] == 1
        q_temp = count_matrix[:, col_idx]

        p2[i, :] = np.nanmean(q_temp, axis=1)

    # Observed co-occurrences
    p3 = np.zeros((num_neurons, num_neurons))
    for i in range(num_neurons):
        for j in range(num_neurons):
            p3[i][j] = np.nanmean(count_matrix[i] * count_matrix[j])
    # remove probability of cell co-occurring with itself
    p3[np.eye(len(p3), dtype=bool)] = np.nan

    # Shuffle
    shuffled_matrix = count_matrix
    num_rows = shuffled_matrix.shape[0]
    num_col = shuffled_matrix.shape[1]

    shuff_p4 = np.zeros((shuffle, num_rows, num_rows))

    for i in range(shuffle):
        this_matrix = shuffled_matrix
        for j in range(num_rows):
            this_matrix[j] = this_matrix[j, np.random.permutation(range(num_col))]
        for k in range(num_rows):
            for m in range(num_rows):
                shuff_p4[i, k, m] = np.nanmean(this_matrix[k]*this_matrix[m])

    # Compare co-occurrences with shuffle
    p4 = np.zeros((num_neurons, num_neurons))
    for i in range(num_neurons):
        for j in range(num_neurons):
            p4[i][j] = (p3[i][j] - np.nanmean(np.squeeze(shuff_p4[:, i, j]))) / np.nanstd(np.squeeze(shuff_p4[:, i, j]))

    p5 = shuff_p4

    # NOT IMPLEMENTED: Handle mask for neurons that were recorded on the same tetrode

    output = dict()
    output['p1'] = p1
    output['p2'] = p2
    output['p3'] = p3
    output['p4'] = p4

    def vector_from_array(array):
        # Old Matlab code indexes by column (aka.Fortran-style), so to get the indices
        # # of the top triangle, we have to do some reshaping.
        # Otherwise, if the vector made up by rows is OK, then simply :
        # # triangle = np.triu_indices(array.size, k=1), out = array[triangle]
        triangle_lower = np.tril_indices(array.shape[0], k=-1)
        flatten_idx = np.arange(array.size).reshape(array.shape)[triangle_lower]
        triangle = np.unravel_index(flatten_idx, array.shape, order='F')

        # triangle = np.triu_indices(array.size, k=1)
        # out = array[triangle]

        return array[triangle]

    for key in output:
        output[key] = vector_from_array(output[key])

    # These probabilities don't need the vector_from_array treatment (p0 is already a vector).
    output['p0'] = p0
    output['p5'] = p5

    return output
