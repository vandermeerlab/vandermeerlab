import vdmlab as vdm

from maze_functions import get_trial_idx, get_zones
from plotting_functions import plot_proportions, plot_bydurations, plot_bytrial

import sys
sys.path.append('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
# sys.path.append('E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
import R063d2_info as r063d2
import R063d3_info as r063d3
import R063d4_info as r063d4
import R063d5_info as r063d5
import R063d6_info as r063d6
import R066d1_info as r066d1
import R066d2_info as r066d2
import R066d4_info as r066d4


# filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\behaviour\\'
filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\behaviour\\'

# infos = [r063d2, r063d3]
infos = [r063d2, r063d3, r063d4, r063d5, r063d6, r066d1, r066d2, r066d4]


durations = dict(u=[], shortcut=[], novel=[])
num_sessions = 0

trials = []

for info in infos:
    print(info.session_id)
    t_start = info.task_times['phase3'][0]
    t_stop = info.task_times['phase3'][1]

    # Slicing events to only Phase 3
    events = info.get_events()
    sliced_events = dict()
    sliced_events['feeder1'] = vdm.time_slice(events['feeder1'], t_start, t_stop)
    sliced_events['feeder2'] = vdm.time_slice(events['feeder2'], t_start, t_stop)

    feeder1_times = sliced_events['feeder1']
    feeder2_times = sliced_events['feeder2']

    spike_pos = get_zones(info, t_start, t_stop)

    trials_idx = get_trial_idx(spike_pos['u']['time'], spike_pos['shortcut']['time'], spike_pos['novel']['time'],
                               feeder1_times, feeder2_times, t_start, t_stop)

    trials.append(trials_idx)

    num_sessions += 1

    for key in durations:
        for trial in trials_idx[key]:
            durations[key].append(trials_idx['stop_trials'][trial[0]] - trials_idx['start_trials'][trial[0]])

durations['num_sessions'] = num_sessions


shortcuts = []
us = []
novels = []
togethers = []

for trial in trials:
    shortcuts.append(len(trial['shortcut'])/float(len(trial['start_trials'])))
    us.append(len(trial['u'])/float(len(trial['start_trials'])))
    novels.append(len(trial['novel'])/float(len(trial['start_trials'])))
    togethers.append(sorted(trial['u'] + trial['shortcut'] + trial['novel']))


plot_proportions(us, shortcuts, novels, filepath)
plot_bydurations(durations, filepath)
plot_bytrial(togethers, filepath)
