from startup import load_csc, load_position, load_videotrack, load_events, load_spikes, convert_to_cm

session_id = 'R042d3'


def get_csc():
    return load_csc('C:\\Users\\Emily\\Desktop\\R042-2013-08-18\\analysis\\cscs\\R042-2013-08-18-csc.mat')


def get_pos(pxl_to_cm):
    pos = load_videotrack('C:\\Users\\Emily\\Desktop\\R042-2013-08-18\\analysis\\positions\\R042-2013-08-18-vt.mat')
    pos['x'] = pos['x'] / pxl_to_cm[0]
    pos['y'] = pos['y'] / pxl_to_cm[1]
    return pos


def get_events():
    return load_events('C:\\Users\\Emily\\Desktop\\R042-2013-08-18\\analysis\\events\\R042-2013-08-18-event.mat')


def get_spikes():
    return load_spikes('C:\\Users\\Emily\\Desktop\\R042-2013-08-18\\analysis\\spikes\\R042-2013-08-18-spike.mat')

# Experimental session-specific task times for R042 day 3 *Note: Alyssa's motivational T experiment
task_times = dict()
task_times['prerecord'] = [2126.64553, 3214.07253]
task_times['task'] = [3238.67853, 5645.16153]
task_times['postrecord'] = [5656.35353, 6563.46453]

pxl_to_cm = (2.9176, 2.3794)

fs = 2000

good_lfp = ['R042-2013-08-18-CSC11a.ncs']

