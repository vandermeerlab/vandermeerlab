import numpy as np

import vdmlab as vdm


def test_swr():
    fs = 2000
    toy_time = np.arange(0, 0.5, 1./fs)
    freq = np.ones(len(toy_time))*100
    freq[int(len(toy_time)*0.4):int(len(toy_time)*0.6)] = 180
    freq[int(len(toy_time)*0.7):int(len(toy_time)*0.9)] = 260
    toy_lfp = np.sin(2.*np.pi*freq*toy_time)

    toy_csc = dict()
    toy_csc['time'] = toy_time
    toy_csc['data'] = toy_lfp

    toy_times, toy_idx, toy_butter = vdm.detect_swr_hilbert(toy_csc, power_thres=0.5, z_thres=0.4)
    assert np.allclose(toy_times['start'][0], 0.1995)
    assert np.allclose(toy_times['stop'][0], 0.3005)
