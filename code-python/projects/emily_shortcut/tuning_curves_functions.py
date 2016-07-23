from shapely.geometry import Point, LineString

import vdmlab as vdm


def linearize(info, pos, t_start, t_stop, expand_by=6):
    # Slicing position to only Phase 3
    t_start_idx = vdm.find_nearest_idx(pos['time'], t_start)
    t_end_idx = vdm.find_nearest_idx(pos['time'], t_stop)
    sliced_pos = dict()
    sliced_pos['x'] = pos['x'][t_start_idx:t_end_idx]
    sliced_pos['y'] = pos['y'][t_start_idx:t_end_idx]
    sliced_pos['time'] = pos['time'][t_start_idx:t_end_idx]

    u_line = LineString(info.u_trajectory)
    shortcut_line = LineString(info.shortcut_trajectory)
    novel_line = LineString(info.novel_trajectory)

    u_start = Point(info.u_trajectory[0])
    u_stop = Point(info.u_trajectory[-1])
    shortcut_start = Point(info.shortcut_trajectory[0])
    shortcut_stop = Point(info.shortcut_trajectory[-1])
    novel_start = Point(info.novel_trajectory[0])
    novel_stop = Point(info.novel_trajectory[-1])
    pedestal_center = Point(info.path_pts['pedestal'][0], info.path_pts['pedestal'][1])
    pedestal = pedestal_center.buffer(expand_by*2.2)

    def expand_line(start_pt, stop_pt, line, expand_by):
        line_expanded = line.buffer(expand_by)
        zone = start_pt.union(line_expanded).union(stop_pt)
        return zone

    zone = dict()
    zone['u'] = expand_line(u_start, u_stop, u_line, expand_by)
    zone['shortcut'] = expand_line(shortcut_start, shortcut_stop, shortcut_line, expand_by)
    zone['novel'] = expand_line(novel_start, novel_stop, novel_line, expand_by)
    zone['ushort'] = zone['u'].intersection(zone['shortcut'])
    zone['unovel'] = zone['u'].intersection(zone['novel'])
    zone['uped'] = zone['u'].intersection(pedestal)
    zone['shortped'] = zone['shortcut'].intersection(pedestal)
    zone['novelped'] = zone['novel'].intersection(pedestal)
    zone['pedestal'] = pedestal

    u_idx = []
    shortcut_idx = []
    novel_idx = []
    other_idx = []
    for pos_idx in range(len(sliced_pos['time'])):
        point = Point([sliced_pos['x'][pos_idx], sliced_pos['y'][pos_idx]])
        if zone['u'].contains(point) or zone['ushort'].contains(point) or zone['unovel'].contains(point):
            u_idx.append(pos_idx)
        elif zone['shortcut'].contains(point) or zone['shortped'].contains(point):
            shortcut_idx.append(pos_idx)
        elif zone['novel'].contains(point) or zone['novelped'].contains(point):
            novel_idx.append(pos_idx)
        else:
            other_idx.append(pos_idx)

    u_pos = vdm.idx_in_pos(sliced_pos, u_idx)
    shortcut_pos = vdm.idx_in_pos(sliced_pos, shortcut_idx)
    novel_pos = vdm.idx_in_pos(sliced_pos, novel_idx)
    other_pos = vdm.idx_in_pos(sliced_pos, other_idx)

    linear = dict()
    linear['u'] = vdm.linear_trajectory(u_pos, u_line, t_start, t_stop)
    linear['shortcut'] = vdm.linear_trajectory(shortcut_pos, shortcut_line, t_start, t_stop)
    linear['novel'] = vdm.linear_trajectory(novel_pos, novel_line, t_start, t_stop)

    return linear, zone
