import numpy as np
from scipy import signal
import scipy.stats as stats
from scipy.signal.signaltools import _next_regular

import matplotlib.pyplot as plt

def butter_bandpass(lowcut, highcut, fs, lfp, order=4):
    """ Filters signal using butterworth filter

    Parameters
    ----------
    lowcut : float
        Suggested 140.0 for sharp-wave ripple detection.
    highcut : float
        Suggested 250.0 for sharp-wave ripple detection.
    fs : int
        Eg. 2000. Should get this from experiment-specifics.
    lfp : np.array
        Eg. csc['data']
    order : int
        Default set to 4.

    Returns
    -------
    filtered_butter : np.array

    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    filtered_butter = signal.filtfilt(b, a, lfp)
    return filtered_butter


def detect_swr_hilbert(csc, fs=2000, lowcut=140.0, highcut=250.0,
                       z_thres=3, power_thres=3, merge_thres=0.02, min_length=0.01):
    """ Finds sharp-wave ripple (SWR) times and indices.

    Parameters
    ----------
    csc : dict
        With time(np.array), data(np.array) as keys
    fs : int
        Experiment-specific, something in the range of 2000 is not unexpected.
    lowcut : float
        The default is set to 140.0
    highcut : float
        The default is set to 250.0
    z_thres : int or float
        The default is set to 5
    power_thres : int or float
        The default is set to 4
    merge_thres : int or float
        The default is set to 0.0
    min_length : float
        Any sequence less than this amount is not considered a sharp-wave ripple.
        The default is set to 0.02.

    Returns
    -------
    swr_times : dict
        With start(float), stop(float) as keys
    swr_idx : dict
        With start(int), stop(int) as keys
    filtered_butter : np.array
        Mostly for plotting purposes

    """
    n_samples = len(csc['data'])

    # Filtering signal with butterworth fitler
    filtered_butter = butter_bandpass(lowcut, highcut, fs, csc['data'])

    # Get LFP power (using Hilbert) and z-score the power
    # Zero padding to nearest regular number to speed up fast fourier transforms (FFT) computed in the hilbert function.
    # Regular numbers are composites of the prime factors 2, 3, and 5.
    hilbert_n = _next_regular(n_samples)
    power_lfp = np.abs(signal.hilbert(filtered_butter, N=hilbert_n))
    power_lfp = power_lfp[:n_samples]  # removing the zero padding now that the power is computed
    zpower_lfp = stats.zscore(power_lfp)

    # Finding locations where the power changes
    detect = zpower_lfp > z_thres
    detect = np.hstack([0, detect, 0])  # pad to detect first or last element change
    signal_change = np.diff(detect.astype(int))

    start_swr_idx = np.where(signal_change == 1)[0]
    stop_swr_idx = np.where(signal_change == -1)[0] - 1

    # Getting times associated with these power changes
    start_time = csc['time'][start_swr_idx]
    stop_time = csc['time'][stop_swr_idx]

    # Merging ranges that are closer - in time - than the merge_threshold.
    no_double = start_time[1:] - stop_time[:-1]
    merge_idx = np.where(no_double < merge_thres)[0]
    start_merged = np.delete(start_time, merge_idx + 1)
    stop_merged = np.delete(stop_time, merge_idx)
    start_merged_idx = np.delete(start_swr_idx, merge_idx + 1)
    stop_merged_idx = np.delete(stop_swr_idx, merge_idx)

    # Removing ranges that are shorter - in time - than the min_length value.
    swr_len = stop_merged - start_merged
    short_idx = np.where(swr_len < min_length)[0]
    start_merged = np.delete(start_merged, short_idx)
    stop_merged = np.delete(stop_merged, short_idx)
    start_merged_idx = np.delete(start_merged_idx, short_idx)
    stop_merged_idx = np.delete(stop_merged_idx, short_idx)

    # Removing ranges that have powers less than the power_threshold if sufficiently different.
    if power_thres > z_thres:
        max_z = []
        for start_idx, stop_idx in zip(start_merged_idx, stop_merged_idx):
            max_z.append(np.max(zpower_lfp[start_idx:stop_idx]))
        max_z = np.array(max_z)

        z_idx = np.where(max_z < power_thres)[0]
        start_merged = np.delete(start_merged, z_idx)
        stop_merged = np.delete(stop_merged, z_idx)
        start_merged_idx = np.delete(start_merged_idx, z_idx)
        stop_merged_idx = np.delete(stop_merged_idx, z_idx)

    swr_idx = dict()
    swr_times = dict()
    swr_times['start'] = start_merged
    swr_times['stop'] = stop_merged
    swr_idx['start'] = start_merged_idx
    swr_idx['stop'] = stop_merged_idx

    print('Number of SWR events found: ', str(len(swr_idx['start'])))

    return swr_times, swr_idx
