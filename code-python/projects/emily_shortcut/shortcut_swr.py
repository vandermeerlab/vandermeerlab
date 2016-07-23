import vdmlab as vdm

from plotting_functions import plot_swrs

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

output_filepath = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\swr\\'
# output_filepath = 'E:\\code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\plots\\swr\\'

for info in infos:
    print('Working on ' + info.session_id)
    csc = info.get_csc()

    swr_times, swr_idx, filtered_butter = vdm.detect_swr_hilbert(csc, fs=info.fs)

    savepath = output_filepath + str(info.session_id) + '-swr_'
    plot_swrs(csc, swr_idx, savepath)
