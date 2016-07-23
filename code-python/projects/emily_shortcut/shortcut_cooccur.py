import os

import vdmlab as vdm

from tuning_curves_functions import get_tc
from plotting_functions import plot_cooccur, plot_solo_coprob

thisdir = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(thisdir, 'info'))
import R063d2_info as r063d2
import R063d3_info as r063d3
import R063d4_info as r063d4
import R063d5_info as r063d5
import R063d6_info as r063d6
import R066d1_info as r066d1
import R066d2_info as r066d2
import R066d4_info as r066d4


pickle_filepath = os.path.join(thisdir, 'cache', 'pickled')
output_filepath = os.path.join(thisdir, 'plots', 'cooccur')

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

        tc = get_tc(info, pos, pickle_filepath)

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

        filename = info.session_id + '_cooccur-' + exp_time + '.png'
        savepath = os.path.join(output_filepath, filename)
        plot_cooccur(probs, savepath)


        # metric = 'p4'
        # title = 'Co-activation above chance levels (p4)'
        # ylabel = 'SWR co-activation z-scored'
        # filename = info.session_id + '_cooccur-' + metric + '_' + exp_time + '.png'
        # savepath = os.path.join(output_filepath, filename)
        # plot_solo_coprob(probs, output_filepath, metric=metric, title=title, ylabel=ylabel, savefig=False)
