from shapely.geometry import Point, LineString

import vdmlab as vdm

import sys
sys.path.append('C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\info')
import R063d2_info as r063d2


output_path = 'C:\\Users\\Emily\\Code\\vandermeerlab\\code-python\\projects\\emily_shortcut\\cache\\matlab\\spike_pos\\'

info = r063d2

pos = info.get_pos(info.pxl_to_cm)
spikes = info.get_spikes()

u_line = LineString(info.u_trajectory)
shortcut_line = LineString(info.shortcut_trajectory)
novel_line = LineString(info.novel_trajectory)

expand_by = 6

u_start = Point(info.u_trajectory[0])
u_stop = Point(info.u_trajectory[-1])
shortcut_start = Point(info.shortcut_trajectory[0])
shortcut_stop = Point(info.shortcut_trajectory[-1])
novel_start = Point(info.novel_trajectory[0])
novel_stop = Point(info.novel_trajectory[-1])
pedestal_center = Point(info.path_pts['pedestal'][0], info.path_pts['pedestal'][1])
pedestal = pedestal_center.buffer(expand_by*2.2)

zone = dict()
zone['u'] = vdm.expand_line(u_start, u_stop, u_line, expand_by)
zone['shortcut'] = vdm.expand_line(shortcut_start, shortcut_stop, shortcut_line, expand_by)
zone['novel'] = vdm.expand_line(novel_start, novel_stop, novel_line, expand_by)
zone['ushort'] = zone['u'].intersection(zone['shortcut'])
zone['unovel'] = zone['u'].intersection(zone['novel'])
zone['uped'] = zone['u'].intersection(pedestal)
zone['shortped'] = zone['shortcut'].intersection(pedestal)
zone['novelped'] = zone['novel'].intersection(pedestal)
zone['pedestal'] = pedestal


spike_position = vdm.spikes_by_position(spikes['time'], zone, pos['time'], pos['x'], pos['y'])


savepath = output_path + 'spike_pos_' + info.session_id + '.mat'
vdm.save_spike_position(spike_position, info.session_id, savepath)
