import pickle
import os

import vdmlab as vdm

from maze_functions import spikes_by_position
from tuning_curves_functions import linearize
from plotting_functions import plot_cooccur

import sys
# sys.path.append('E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
sys.path.append('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
import R063d2_info as r063d2
import R063d3_info as r063d3
import R063d4_info as r063d4
import R063d5_info as r063d5
import R063d6_info as r063d6
import R066d1_info as r066d1
import R066d2_info as r066d2
import R066d4_info as r066d4


# pickle_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\pickled\\'
# output_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\output\\'
pickle_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\pickled\\'
output_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\cooccur\\'


# info = R063d3
infos = [r063d2, r063d3, r063d4, r063d5, r063d6, r066d1, r066d2, r066d4]

exp_times = ['pauseA', 'pauseB']
for info in infos:
    for exp_time in exp_times:
        csc = info.get_csc()
        pos = info.get_pos(info.pxl_to_cm)
        spikes = info.get_spikes()

        t_start = info.task_times[exp_time][0]
        t_stop = info.task_times[exp_time][1]

        t_start_idx = vdm.find_nearest_idx(csc['time'], t_start)
        t_end_idx = vdm.find_nearest_idx(csc['time'], t_stop)

        sliced_csc = dict()
        sliced_csc['time'] = csc['time'][t_start_idx:t_end_idx]
        sliced_csc['data'] = csc['data'][t_start_idx:t_end_idx]

        swr_times, swr_idx, filtered_butter = vdm.detect_swr_hilbert(sliced_csc, fs=info.fs, power_thres=5)

        # pickled_spike_pos = pickle_filepath + info.session_id + '_spike_position_phase3.pkl'
        # if os.path.isfile(pickled_spike_pos):
        #     with open(pickled_spike_pos, 'rb') as fileobj:
        #         spike_position = pickle.load(fileobj)
        # else:
        #     t_start = info.task_times['phase3'][0]
        #     t_stop = info.task_times['phase3'][1]
        #
        #     spikes = info.get_spikes()
        #     pos = info.get_pos(info.pxl_to_cm)
        #
        #     linear, zone = linearize(info, pos, t_start, t_stop)
        #     spike_position = spikes_by_position(spikes['time'], zone, pos['time'], pos['x'], pos['y'])
        #     with open(pickled_spike_pos, 'wb') as fileobj:
        #         pickle.dump(spike_position, fileobj)

        pickled_tc = pickle_filepath + info.session_id + '_tuning_curves_phase3.pkl'
        if os.path.isfile(pickled_tc):
            with open(pickled_tc, 'rb') as fileobj:
                tc = pickle.load(fileobj)
        else:
            t_start = info.task_times['phase3'][0]
            t_stop = info.task_times['phase3'][1]

            # spikes = info.get_spikes()
            pos = info.get_pos(info.pxl_to_cm)

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

        all_u_fields = vdm.find_fields(tc['u'])
        all_shortcut_fields = vdm.find_fields(tc['shortcut'])
        all_novel_fields = vdm.find_fields(tc['novel'])

        u_fields_unique = vdm.unique_fields(all_u_fields, all_shortcut_fields, all_novel_fields)
        shortcut_fields_unique = vdm.unique_fields(all_shortcut_fields, all_u_fields, all_novel_fields)
        novel_fields_unique = vdm.unique_fields(all_novel_fields, all_u_fields, all_shortcut_fields)

        u_fields_single = vdm.get_single_field(u_fields_unique)
        shortcut_fields_single = vdm.get_single_field(shortcut_fields_unique)
        novel_fields_single = vdm.get_single_field(novel_fields_unique)

        u_spikes = []
        for key in u_fields_unique:
            u_spikes.append(spikes['time'][key])

        shortcut_spikes = []
        for key in shortcut_fields_unique:
            shortcut_spikes.append(spikes['time'][key])

        novel_spikes = []
        for key in novel_fields_unique:
            novel_spikes.append(spikes['time'][key])

        count_matrix = dict()
        count_matrix['u'] = vdm.spike_counts(u_spikes, swr_times, window=0.1)
        count_matrix['shortcut'] = vdm.spike_counts(shortcut_spikes, swr_times, window=0.1)
        count_matrix['novel'] = vdm.spike_counts(novel_spikes, swr_times, window=0.1)

        probs = dict()
        probs['u'] = vdm.cooccur_probabilities(count_matrix['u'])
        probs['shortcut'] = vdm.cooccur_probabilities(count_matrix['shortcut'])
        probs['novel'] = vdm.cooccur_probabilities(count_matrix['novel'])

        # metric_title = OrderedDict()
        # metric_title['p0'] = ['Probability that a given neuron is active (p0)', 'Proportion of SWRs active']
        # metric_title['p3'] = ['Observed co-activity (p3)', 'Cell pair joint probability']
        # metric_title['p4'] = ['Co-activation above chance levels (p4)', 'SWR co-activation z-scored']

        cooccur_filename = output_filepath + info.session_id + '_cooccur-' + exp_time + '.png'
        plot_cooccur(probs, cooccur_filename)
