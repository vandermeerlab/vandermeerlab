import os

from startup import load_csc, load_videotrack, load_events, load_spikes, convert_to_cm

session_id = 'R063d4'

thisdir = os.path.dirname(os.path.realpath(__file__))
dataloc = os.path.abspath(os.path.join(thisdir, '..', 'cache', 'data'))

species = 'rat'
behavior = 'shortcut'
target = 'dCA1'
experimenter = 'Emily Irvine'


def get_csc():
    return load_csc(os.path.join(dataloc, 'R063-2015-03-23-csc.mat'))


def get_pos(pxl_to_cm):
    pos = load_videotrack(os.path.join(dataloc, 'R063-2015-03-23-vt.mat'))
    pos['x'] = pos['x'] / pxl_to_cm[0]
    pos['y'] = pos['y'] / pxl_to_cm[1]
    return pos


def get_events():
    return load_events(os.path.join(dataloc, 'R063-2015-03-23-event.mat'))


def get_spikes():
    return load_spikes(os.path.join(dataloc, 'R063-2015-03-23-spike.mat'))


# Experimental session-specific task times for R063 day 4
task_times = dict()
task_times['prerecord'] = [1074.6, 1378.7]
task_times['phase1'] = [1415.9, 1847.2]
task_times['pauseA'] = [1860.6, 2486.0]
task_times['phase2'] = [2504.6, 3704.5]
task_times['pauseB'] = [3725.3, 5600.7]
task_times['phase3'] = [5627.4, 8638.8]
task_times['postrecord'] = [8656.4, 9000.0]

pxl_to_cm = (7.9628, 7.2755)

fs = 2000

good_lfp = ['R063-2015-03-23-CSC11b.ncs']
good_swr = ['']
good_theta = ['']

# Session-specific path trajectory points
path_pts = dict()
path_pts['feeder1'] = [552, 461]
path_pts['point1'] = [551, 409]
path_pts['point2'] = [547, 388]
path_pts['point3'] = [535, 383]
path_pts['point4'] = [479, 381]
path_pts['point5'] = [370, 400]
path_pts['point6'] = [249, 371]
path_pts['point7'] = [206, 378]
path_pts['point8'] = [219, 325]
path_pts['point9'] = [234, 98]
path_pts['point10'] = [244, 66]
path_pts['point11'] = [275, 51]
path_pts['point12'] = [334, 46]
path_pts['feeder2'] = [662, 55]
path_pts['shortcut1'] = [546, 387]
path_pts['point13'] = [620, 373]
path_pts['point14'] = [648, 390]
path_pts['point15'] = [659, 355]
path_pts['point16'] = [661, 150]
path_pts['shortcut2'] = [662, 55]
path_pts['novel1'] = [331, 392]
path_pts['point17'] = [324, 460]
path_pts['point18'] = [316, 471]
path_pts['point19'] = [289, 478]
path_pts['novel2'] = [113, 470]
path_pts['pedestal'] = [412, 203]

path_pts = convert_to_cm(path_pts, pxl_to_cm)

u_trajectory = [path_pts['feeder1'], path_pts['point1'], path_pts['point2'],
                path_pts['point3'], path_pts['point4'], path_pts['point5'],
                path_pts['point6'], path_pts['point7'], path_pts['point8'],
                path_pts['point9'], path_pts['point10'], path_pts['point11'],
                path_pts['point12'], path_pts['feeder2']]

shortcut_trajectory = [path_pts['shortcut1'], path_pts['point13'], path_pts['point14'],
                       path_pts['point15'], path_pts['point16'], path_pts['shortcut2']]

novel_trajectory = [path_pts['novel1'], path_pts['point17'], path_pts['point18'],
                    path_pts['point19'], path_pts['novel2']]

sequence = dict()
sequence['swr_start'] = [7946.591335, 5101.158335]
sequence['swr_stop'] = [7946.741335, 5101.233835]
sequence['run_start'] = [5817.4, 2882.6]
sequence['run_stop'] = [5837.4, 2912.6]
sequence['ms'] = 10
sequence['loc'] = 1
sequence['colours'] = ['#bd0026', '#fc4e2a', '#ef3b2c', '#ec7014', '#fe9929', '#78c679',
                       '#41ab5d', '#238443', '#66c2a4', '#41b6c4', '#1d91c0', '#8c6bb1']
