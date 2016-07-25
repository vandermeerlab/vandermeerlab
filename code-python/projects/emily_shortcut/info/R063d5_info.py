import os

from startup import load_csc, load_videotrack, load_events, load_spikes, convert_to_cm

session_id = 'R063d5'

thisdir = os.path.dirname(os.path.realpath(__file__))
dataloc = os.path.abspath(os.path.join(thisdir, '..', 'cache', 'data'))

species = 'rat'
behavior = 'shortcut'
target = 'dCA1'
experimenter = 'Emily Irvine'


def get_csc():
    return load_csc(os.path.join(dataloc, 'R063-2015-03-24-csc.mat'))


def get_pos(pxl_to_cm):
    pos = load_videotrack(os.path.join(dataloc, 'R063-2015-03-24-vt.mat'))
    pos['x'] = pos['x'] / pxl_to_cm[0]
    pos['y'] = pos['y'] / pxl_to_cm[1]
    return pos


def get_events():
    return load_events(os.path.join(dataloc, 'R063-2015-03-24-event.mat'))


def get_spikes():
    return load_spikes(os.path.join(dataloc, 'R063-2015-03-24-spike.mat'))

# Experimental session-specific task times for R063 day 5
task_times = dict()
task_times['prerecord'] = [2160.2, 2479.9]
task_times['phase1'] = [2516.0, 3119.6]
task_times['pauseA'] = [3136.2, 3768.3]
task_times['phase2'] = [3797.6, 5066.1]
task_times['pauseB'] = [5088.4, 7189.4]
task_times['phase3'] = [7221.2, 10264.0]
task_times['postrecord'] = [10285.0, 10591.0]

pxl_to_cm = (7.4906, 7.2379)

fs = 2000

good_lfp = ['R063-2015-03-24-CSC13d.ncs']
good_swr = ['']
good_theta = ['']

# Session-specific path trajectory points
path_pts = dict()
path_pts['feeder1'] = [551, 467]
path_pts['point1'] = [549, 413]
path_pts['point2'] = [538, 380]
path_pts['point3'] = [507, 379]
path_pts['point4'] = [390, 399]
path_pts['point5'] = [258, 375]
path_pts['point6'] = [226, 368]
path_pts['point7'] = [217, 330]
path_pts['point8'] = [220, 235]
path_pts['point9'] = [232, 97]
path_pts['point10'] = [240, 67]
path_pts['point11'] = [277, 51]
path_pts['point12'] = [517, 42]
path_pts['feeder2'] = [658, 66]

path_pts['shortcut1'] = [328, 384]
path_pts['point13'] = [333, 191]
path_pts['point14'] = [350, 158]
path_pts['point15'] = [386, 144]
path_pts['point16'] = [419, 138]
path_pts['point17'] = [441, 141]
path_pts['point18'] = [447, 100]
path_pts['shortcut2'] = [451, 44]

path_pts['novel1'] = [234, 80]
path_pts['point19'] = [152, 79]
path_pts['point20'] = [130, 93]
path_pts['point21'] = [127, 117]
path_pts['novel2'] = [124, 272]
path_pts['pedestal'] = [562, 233]

path_pts = convert_to_cm(path_pts, pxl_to_cm)

u_trajectory = [path_pts['feeder1'], path_pts['point1'], path_pts['point2'],
                path_pts['point3'], path_pts['point4'], path_pts['point5'],
                path_pts['point6'], path_pts['point7'], path_pts['point8'],
                path_pts['point9'], path_pts['point10'], path_pts['point11'],
                path_pts['point12'], path_pts['feeder2']]

shortcut_trajectory = [path_pts['shortcut1'], path_pts['point13'], path_pts['point14'],
                       path_pts['point15'], path_pts['point16'], path_pts['point17'],
                       path_pts['point18'], path_pts['shortcut2']]

novel_trajectory = [path_pts['novel1'], path_pts['point19'], path_pts['point20'],
                    path_pts['point21'], path_pts['novel2']]

sequence = dict()
sequence['swr_start'] = [5730.325233, 6376.529733, 6967.243733]
sequence['swr_stop'] = [5730.413733, 6376.624233, 6967.307233]
sequence['run_start'] = [3997.6, 4062.6, 3892.6]
sequence['run_stop'] = [4027.6, 4092.6, 3922.6]
sequence['ms'] = 10
sequence['loc'] = 1
sequence['colours'] = ['#bd0026', '#fc4e2a', '#ef3b2c', '#ec7014', '#fe9929', '#78c679',
                       '#41ab5d', '#238443', '#66c2a4', '#41b6c4', '#1d91c0']
