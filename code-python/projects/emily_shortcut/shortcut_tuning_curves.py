import pickle
import os

import vdmlab as vdm

from maze_functions import spikes_by_position
from tuning_curves_functions import linearize

import sys
sys.path.append('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
# sys.path.append('E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
import R042d3_info as r042d3
import R063d2_info as r063d2
import R063d3_info as r063d3
import R063d4_info as r063d4
import R063d5_info as r063d5
import R063d6_info as r063d6
import R066d1_info as r066d1
import R066d2_info as r066d2
import R066d4_info as r066d4


# infos = [r042d3]
# infos = [r063d2]
infos = [r063d2, r063d3, r063d4, r063d5, r063d6, r066d1, r066d2, r066d4]

# pickle_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\pickled\\'
# output_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\'
pickle_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\pickled\\'
output_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\'

for info in infos:
    pos = info.get_pos(info.pxl_to_cm)

    pickled_tc = pickle_filepath + info.session_id + '_tuning_curves_phase3.pkl'
    if os.path.isfile(pickled_tc):
        with open(pickled_tc, 'rb') as fileobj:
            tc = pickle.load(fileobj)
    else:
        t_start = info.task_times['phase3'][0]
        t_stop = info.task_times['phase3'][1]

        spikes = info.get_spikes()

        linear, zone = linearize(info, pos, t_start, t_stop)

        pickled_spike_pos = pickle_filepath + info.session_id + '_spike_position_phase3.pkl'
        if os.path.isfile(pickled_spike_pos):
            with open(pickled_spike_pos, 'rb') as fileobj:
                spike_position = pickle.load(fileobj)
        else:
            sliced_spikes = vdm.time_slice(spikes['time'], t_start, t_stop)
            spike_position = spikes_by_position(sliced_spikes, zone, pos['time'], pos['x'], pos['y'])
            with open(pickled_spike_pos, 'wb') as fileobj:
                pickle.dump(spike_position, fileobj)

        tc = dict()
        tc['u'] = vdm.tuning_curve(linear['u'], spike_position['u'], num_bins=47)
        tc['shortcut'] = vdm.tuning_curve(linear['shortcut'], spike_position['shortcut'], num_bins=47)
        tc['novel'] = vdm.tuning_curve(linear['novel'], spike_position['novel'], num_bins=47)
        with open(pickled_tc, 'wb') as fileobj:
            pickle.dump(tc, fileobj)
