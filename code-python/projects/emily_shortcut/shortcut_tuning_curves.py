import os
import vdmlab as vdm

from tuning_curves_functions import get_tc, get_odd_firing_idx
from plotting_functions import plot_sorted_tc

import info.R063d2_info as r063d2
import info.R063d3_info as r063d3
import info.R063d4_info as r063d4
import info.R063d5_info as r063d5
import info.R063d6_info as r063d6
import info.R066d1_info as r066d1
import info.R066d2_info as r066d2
import info.R066d4_info as r066d4


thisdir = os.path.dirname(os.path.realpath(__file__))


# infos = [r063d2]
infos = [r063d2, r063d3, r063d4, r063d5, r063d6, r066d1, r066d2, r066d4]

pickle_filepath = os.path.join(thisdir, 'cache', 'pickled')
output_filepath = os.path.join(thisdir, 'plots', 'tuning')

for info in infos:
    print(info.session_id)
    pos = info.get_pos(info.pxl_to_cm)

    tc = get_tc(info, pos, pickle_filepath)

    sort_idx = dict()
    odd_firing_idx = dict()
    sorted_tc = dict(u=[], shortcut=[], novel=[])

    for key in tc:
        sort_idx[key] = vdm.get_sort_idx(tc[key])
        odd_firing_idx[key] = get_odd_firing_idx(tc[key])

        for idx in sort_idx[key]:
            if idx not in odd_firing_idx[key]:
                sorted_tc[key].append(tc[key][idx])

        filename = info.session_id + '-sorted_tc-' + key + '.png'
        savepath = os.path.join(output_filepath, filename)
        plot_sorted_tc(sorted_tc[key], savepath)
