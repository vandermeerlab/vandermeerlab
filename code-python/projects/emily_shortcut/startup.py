import scipy.io as sio
import numpy as np

# Loading functions from *.mat files
def load_csc(matfile):
    loading_csc = sio.loadmat(matfile)
    csc = dict(time=[])
    csc['data'] = loading_csc['csc_data'][0]
    csc['time'] = loading_csc['csc_tvec']
    csc['type'] = loading_csc['csc_type'][0]
    csc['label'] = loading_csc['csc_label'][0][0][0]
    return csc


def load_position(matfile):
    loading_pos = sio.loadmat(matfile)
    pos = dict()
    pos['x'] = loading_pos['pos_datax'][0]
    pos['y'] = loading_pos['pos_datay'][0]
    pos['time'] = loading_pos['pos_tvec'][0]
    pos['type'] = loading_pos['pos_type'][0]
    pos['label'] = loading_pos['pos_label'][0][0][0]
    return pos

# This data had issues with the feeder lights contaminating the position tracking, so those contaminating
# signals were removed.


def load_videotrack(matfile):
    loading_vt = sio.loadmat(matfile)
    vt = dict()
    vt['time'] = loading_vt['pos_tsd'][0][0][1][0]
    # vt['x'] = loading_vt['pos_tsd'][0][0][1][0]
    # vt['y'] = loading_vt['pos_tsd'][0][0][1][1]
    vt['x'] = loading_vt['pos_tsd'][0][0][2][0]
    vt['y'] = loading_vt['pos_tsd'][0][0][2][1]

    nan_idx = np.isnan(vt['x']) | np.isnan(vt['y'])
    vt['time'] = vt['time'][~nan_idx]
    vt['x'] = vt['x'][~nan_idx]
    vt['y'] = vt['y'][~nan_idx]
    return vt


def load_events(matfile):
    loading_events = sio.loadmat(matfile)
    events = dict()
    events['led1'] = loading_events['evt_led1id'][0]
    events['led2'] = loading_events['evt_led2id'][0]
    events['ledoff'] = loading_events['evt_ledoff'][0]
    events['pb1'] = loading_events['evt_pb1id'][0]
    events['pb2'] = loading_events['evt_pb2id'][0]
    events['pboff'] = loading_events['evt_pboff'][0]
    events['feeder1'] = loading_events['evt_feeder1id'][0]
    events['feeder2'] = loading_events['evt_feeder2id'][0]
    events['feederoff'] = loading_events['evt_feederoff'][0]
    events['type'] = loading_events['evt_type'][0]
    events['label'] = loading_events['evt_label'][0]
    return events


def load_spikes(matfile):
    loading_spikes = sio.loadmat(matfile)
    spikes = dict()
    spikes['time'] = loading_spikes['spikes_times'][0]
    spikes['type'] = loading_spikes['spikes_type'][0]
    spikes['label'] = loading_spikes['spikes_label'][0]
    return spikes


def convert_to_cm(path_pts, xy_conversion):
    for key in path_pts:
        path_pts[key][0] = path_pts[key][0] / xy_conversion[0]
        path_pts[key][1] = path_pts[key][1] / xy_conversion[1]
    return path_pts