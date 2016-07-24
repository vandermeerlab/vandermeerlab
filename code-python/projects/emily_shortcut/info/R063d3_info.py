import os

from startup import load_csc, load_position, load_videotrack, load_events, load_spikes, convert_to_cm

session_id = 'R063d3'
dataloc = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\data\\'

species = 'rat'
behavior = 'shortcut'
target = 'dCA1'
experimenter = 'Emily Irvine'

def get_csc():
    return load_csc(os.path.join(dataloc, 'R063-2015-03-22-csc.mat'))


def get_pos(pxl_to_cm):
    pos = load_videotrack(os.path.join(dataloc, 'R063-2015-03-22-vt.mat'))
    pos['x'] = pos['x'] / pxl_to_cm[0]
    pos['y'] = pos['y'] / pxl_to_cm[1]
    return pos


def get_events():
    return load_events(os.path.join(dataloc, 'R063-2015-03-22-event.mat'))


def get_spikes():
    return load_spikes(os.path.join(dataloc, 'R063-2015-03-22-spike.mat'))

# plt.plot(pos['x'], pos['y'])
# plt.show()

# Experimental session-specific task times for R063 day 3
task_times = dict()
task_times['prerecord'] = [837.4714, 1143.1]
task_times['phase1'] = [1207.9, 2087.5]
task_times['pauseA'] = [2174.3, 2800.8]
task_times['phase2'] = [2836.2, 4034.1]
task_times['pauseB'] = [4051.3, 6185.6]
task_times['phase3'] = [6249.5, 9373.7]
task_times['postrecord'] = [9395.4, 9792.5]

pxl_to_cm = (7.3452, 7.2286)

fs = 2000

good_lfp = ['R063-2015-03-22-CSC14b.ncs']
good_swr = ['']
good_theta = ['']

# Session-specific path trajectory points
path_pts = dict()
path_pts['point1'] = [567, 463]
path_pts['feeder1'] = [547, 469]
path_pts['point2'] = [542, 396]
path_pts['turn1'] = [536, 377]
path_pts['point3'] = [511, 380]
path_pts['point4'] = [453, 408]
path_pts['point5'] = [316, 400]
path_pts['point6'] = [335, 425]
path_pts['point7'] = [389, 402]
path_pts['point8'] = [248, 370]
path_pts['turn2'] = [217, 375]
path_pts['point9'] = [217, 316]
path_pts['point10'] = [236, 84]
path_pts['turn3'] = [249, 59]
path_pts['point11'] = [289, 51]
path_pts['point12'] = [532, 47]
path_pts['feeder2'] = [670, 56]
path_pts['shortcut1'] = [446, 391]
path_pts['point13'] = [438, 334]
path_pts['point14'] = [449, 295]
path_pts['point15'] = [471, 277]
path_pts['point16'] = [621, 269]
path_pts['point17'] = [649, 263]
path_pts['point18'] = [666, 280]
path_pts['point19'] = [660, 240]
path_pts['shortcut2'] = [654, 56]
path_pts['novel1'] = [247, 61]
path_pts['point20'] = [135, 55]
path_pts['point21'] = [128, 64]
path_pts['point22'] = [132, 83]
path_pts['novel2'] = [130, 266]
path_pts['pedestal'] = [371, 172]

path_pts = convert_to_cm(path_pts, pxl_to_cm)

u_trajectory = [path_pts['point1'], path_pts['feeder1'], path_pts['point2'], path_pts['turn1'],
                path_pts['point3'], path_pts['point4'], path_pts['point5'], path_pts['point6'],
                path_pts['point7'], path_pts['point8'], path_pts['turn2'], path_pts['point9'],
                path_pts['point10'], path_pts['turn3'], path_pts['point11'], path_pts['point12'],
                path_pts['feeder2']]

shortcut_trajectory = [path_pts['shortcut1'], path_pts['point13'], path_pts['point14'],
                       path_pts['point15'], path_pts['point16'], path_pts['point17'], path_pts['point18'],
                       path_pts['point19'], path_pts['shortcut2']]

novel_trajectory = [path_pts['novel1'], path_pts['point20'], path_pts['point21'],
                    path_pts['point22'], path_pts['novel2']]


