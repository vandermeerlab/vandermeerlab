from startup import load_csc, load_position, load_videotrack, load_events, load_spikes, convert_to_cm

session_id = 'R066d2'


def get_csc():
    return load_csc('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\'
                    'cscs\\R066-2014-11-28-csc.mat')


def get_pos(pxl_to_cm):
    pos = load_videotrack('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\'
                          'positions\\R066-2014-11-28-vt.mat')
    pos['x'] = pos['x'] / pxl_to_cm[0]
    pos['y'] = pos['y'] / pxl_to_cm[1]
    return pos


def get_events():
    return load_events('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\'
                       'events\\R066-2014-11-28-event.mat')


def get_spikes():
    return load_spikes('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\'
                       'spikes\\R066-2014-11-28-spike.mat')


# plt.plot(pos['x'], pos['y'], 'b.', ms=1)
# plt.show()

# Experimental session-specific task times for R066 day 2
task_times = dict()
task_times['prerecord'] = [11850.0, 12155.0]
task_times['phase1'] = [12210.0, 12840.0]
task_times['pauseA'] = [12900.0, 13501.0]
task_times['phase2'] = [13574.0, 14776.0]
task_times['pauseB'] = [14825.0, 16633.0]
task_times['phase3'] = [16684.0, 19398.0]
task_times['postrecord'] = [19436.0, 19742.0]

pxl_to_cm = (7.5460, 7.2192)

fs = 2000

good_lfp = ['R066-2014-11-28-CSC11d.ncs']

# Session-specific path trajectory points
path_pts = dict()
path_pts['feeder1'] = [530, 460]
path_pts['pt1'] = [525, 382]
path_pts['pt2'] = [472, 375]
path_pts['pt3'] = [425, 397]
path_pts['pt4'] = [404, 359]
path_pts['pt5'] = [396, 396]
path_pts['pt6'] = [348, 395]
path_pts['pt7'] = [307, 357]
path_pts['pt8'] = [298, 390]
path_pts['pt9'] = [266, 370]
path_pts['pt10'] = [222, 360]
path_pts['pt11'] = [204, 365]
path_pts['pt12'] = [194, 304]
path_pts['pt13'] = [207, 88]
path_pts['pt14'] = [206, 51]
path_pts['pt15'] = [269, 44]
path_pts['pt16'] = [536, 48]
path_pts['feeder2'] = [638, 51]
path_pts['pt17'] = [665, 51]
path_pts['shortcut1'] = [525, 382]
path_pts['spt1'] = [530, 203]
path_pts['spt2'] = [550, 173]
path_pts['spt3'] = [532, 168]
path_pts['spt4'] = [630, 178]
path_pts['shortcut2'] = [638, 51]
path_pts['novel1'] = [204, 365]
path_pts['npt1'] = [89, 359]
path_pts['novel2'] = [98, 149]
path_pts['pedestal'] = [331, 206]

path_pts = convert_to_cm(path_pts, pxl_to_cm)

u_trajectory = [path_pts['feeder1'], path_pts['pt1'], path_pts['pt2'],
                path_pts['pt3'], path_pts['pt4'], path_pts['pt5'],
                path_pts['pt6'], path_pts['pt7'], path_pts['pt8'],
                path_pts['pt9'], path_pts['pt10'], path_pts['pt11'],
                path_pts['pt12'], path_pts['pt13'], path_pts['pt14'],
                path_pts['pt15'], path_pts['pt16'], path_pts['feeder2'], path_pts['pt17']]

shortcut_trajectory = [path_pts['shortcut1'], path_pts['spt1'], path_pts['spt2'],
                       path_pts['spt3'], path_pts['spt4'], path_pts['shortcut2']]

novel_trajectory = [path_pts['novel1'], path_pts['npt1'], path_pts['novel2']]


