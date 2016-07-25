import os

from startup import load_csc, load_videotrack, load_events, load_spikes, convert_to_cm

session_id = 'R066d4'

thisdir = os.path.dirname(os.path.realpath(__file__))
dataloc = os.path.abspath(os.path.join(thisdir, '..', 'cache', 'data'))

species = 'rat'
behavior = 'shortcut'
target = 'dCA1'
experimenter = 'Emily Irvine'


def get_csc():
    return load_csc(os.path.join(dataloc, 'R066-2014-12-01-csc.mat'))


def get_pos(pxl_to_cm):
    pos = load_videotrack(os.path.join(dataloc, 'R066-2014-12-01-vt.mat'))
    pos['x'] = pos['x'] / pxl_to_cm[0]
    pos['y'] = pos['y'] / pxl_to_cm[1]
    return pos


def get_events():
    return load_events(os.path.join(dataloc, 'R066-2014-12-01-event.mat'))


def get_spikes():
    return load_spikes(os.path.join(dataloc, 'R066-2014-12-01-spike.mat'))


# Experimental session-specific task times for R066 day 4
task_times = dict()
task_times['prerecord'] = [8.8210e+03, 9.1346e+03]
task_times['phase1'] = [9.1677e+03, 9.6490e+03]
task_times['pauseA'] = [9.7725e+03, 1.0374e+04]
task_times['phase2'] = [1.0406e+04, 1.1606e+04]
task_times['pauseB'] = [1.1675e+04, 1.3479e+04]
task_times['phase3'] = [1.3514e+04, 1.5619e+04]
task_times['postrecord'] = [1.5650e+04, 1.6257e+04]

pxl_to_cm = (7.6032, 7.1722)

fs = 2000

good_lfp = ['R066-2014-12-01-CSC02b.ncs']
good_swr = ['']
good_theta = ['']

# Session-specific path trajectory points
path_pts = dict()
path_pts['feeder1'] = [525, 453]
path_pts['point1a'] = [543, 426]
path_pts['point1'] = [533, 400]
path_pts['point2'] = [539, 380]
path_pts['point3'] = [539, 349]
path_pts['point4'] = [491, 372]
path_pts['point5'] = [436, 368]
path_pts['point6'] = [351, 402]
path_pts['point7'] = [314, 371]
path_pts['point8'] = [278, 402]
path_pts['point9'] = [266, 372]
path_pts['point10'] = [218, 374]
path_pts['point11'] = [194, 378]
path_pts['point12'] = [207, 308]
path_pts['point13'] = [197, 82]
path_pts['point14'] = [210, 54]
path_pts['point15'] = [286, 46]
path_pts['point16'] = [520, 44]
path_pts['feeder2'] = [637, 63]
path_pts['shortcut1'] = [566, 373]
path_pts['point17'] = [602, 368]
path_pts['point18'] = [627, 366]
path_pts['point19'] = [636, 330]
path_pts['point20'] = [630, 134]
path_pts['shortcut2'] = [637, 63]
path_pts['novel1'] = [302, 405]
path_pts['point21'] = [303, 455]
path_pts['point22'] = [316, 471]
path_pts['point23'] = [289, 478]
path_pts['novel2'] = [113, 470]
path_pts['pedestal'] = [337, 225]

path_pts = convert_to_cm(path_pts, pxl_to_cm)

u_trajectory = [path_pts['feeder1'], path_pts['point1'], path_pts['point1a'], path_pts['point2'],
                path_pts['point3'], path_pts['point4'], path_pts['point5'],
                path_pts['point6'], path_pts['point7'], path_pts['point8'],
                path_pts['point9'], path_pts['point10'], path_pts['point11'],
                path_pts['point12'], path_pts['point13'], path_pts['point14'],
                path_pts['point15'], path_pts['point16'], path_pts['feeder2']]

shortcut_trajectory = [path_pts['shortcut1'], path_pts['point17'], path_pts['point18'],
                       path_pts['point19'], path_pts['point20'], path_pts['shortcut2']]

novel_trajectory = [path_pts['novel1'], path_pts['point21'], path_pts['point22'],
                    path_pts['point23'], path_pts['novel2']]


