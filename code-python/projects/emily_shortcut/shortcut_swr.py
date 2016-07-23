import os

import vdmlab as vdm

from plotting_functions import plot_swrs

thisdir = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(os.path.join(thisdir, 'info')))
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
infos = [r063d2]
# infos = [r063d2, r063d3, r063d4, r063d5, r063d6, r066d1, r066d2, r066d4]

output_filepath = os.path.join(thisdir, 'plots', 'swr')


for info in infos:
    print('Working on ' + info.session_id)
    csc = info.get_csc()

    swr_times, swr_idx, filtered_butter = vdm.detect_swr_hilbert(csc, fs=info.fs)

    filename = info.session_id + '-swr_'
    saveloc = os.path.join(output_filepath, filename)
    plot_swrs(csc, swr_idx, saveloc)
