import numpy as np
import os.path
import pickle

import vdmlab as vdm

from maze_functions import spikes_by_position
from tuning_curves_functions import linearize
from plotting_functions import plot_fields

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


# infos = [r066d1]
infos = [r063d2, r063d3, r063d4, r063d5, r063d6, r066d1, r066d2, r066d4]

# pickle_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\pickled\\'
# output_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\fields\\'
pickle_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\pickled\\'
output_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\fields\\'


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


    pickled_spike_heatmaps = pickle_filepath + info.session_id + '_spike_heatmaps.pkl'
    if os.path.isfile(pickled_spike_heatmaps):
        with open(pickle_filepath + info.session_id + '_spike_heatmaps.pkl', 'rb') as fileobj:
            spike_heatmaps = pickle.load(fileobj)
    else:
        spikes = info.get_spikes()

        all_neurons = list(range(0, len(spikes['time'])))
        spike_heatmaps = vdm.get_heatmaps(all_neurons, spikes, pos)
        with open(pickled_spike_heatmaps, 'wb') as fileobj:
            pickle.dump(spike_heatmaps, fileobj)

    all_u_fields = vdm.find_fields(tc['u'])
    all_shortcut_fields = vdm.find_fields(tc['shortcut'])
    all_novel_fields = vdm.find_fields(tc['novel'])

    u_compare = vdm.find_fields(tc['u'], hz_thres=3)
    shortcut_compare = vdm.find_fields(tc['shortcut'], hz_thres=3)
    novel_compare = vdm.find_fields(tc['novel'], hz_thres=3)

    u_fields_unique = vdm.unique_fields(all_u_fields, shortcut_compare, novel_compare)
    shortcut_fields_unique = vdm.unique_fields(all_shortcut_fields, u_compare, novel_compare)
    novel_fields_unique = vdm.unique_fields(all_novel_fields, u_compare, shortcut_compare)

    u_fields = vdm.sized_fields(u_fields_unique)
    shortcut_fields = vdm.sized_fields(shortcut_fields_unique)
    novel_fields = vdm.sized_fields(novel_fields_unique)

    u_fields_single = vdm.get_single_field(u_fields)
    shortcut_fields_single = vdm.get_single_field(shortcut_fields)
    novel_fields_single = vdm.get_single_field(novel_fields)

    print('U: Of', str(len(all_u_fields)), 'fields,',
          str(len(u_fields)), 'are unique, with',
          str(len(u_fields_single)), 'with single peaks.')
    print('Shortcut: Of', str(len(all_shortcut_fields)), 'fields,',
          str(len(shortcut_fields)), 'are unique, with',
          str(len(shortcut_fields_single)), 'with single peaks.')
    print('Novel: Of', str(len(all_novel_fields)), 'fields,',
          str(len(novel_fields)), 'are unique, with',
          str(len(novel_fields_single)), 'with single peaks.')

    num_bins = 100

    all_trajectories = dict(u=u_fields, shortcut=shortcut_fields, novel=novel_fields)

    for trajectory in all_trajectories:
        all_heatmaps = np.zeros((num_bins, num_bins))
        for key in all_trajectories[trajectory]:
            all_heatmaps += spike_heatmaps[key]
        savepath = output_filepath + info.session_id + '-fields_' + str(trajectory) + '.png'
        plot_fields(all_heatmaps, pos, savepath, num=len(all_trajectories[trajectory]))

# for key in novel_fields_unique:
#     print('plotting neuron ' + str(key))
#     plot_fields(spike_heatmaps[key], pos)
