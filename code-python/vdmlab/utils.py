import numpy as np
from matplotlib.offsetbox import AnchoredOffsetbox


def find_nearest_idx(array, val):
    """Finds nearest index in array to value.

    Parameters
    ----------
    array : numpy array
    val : float

    Returns
    -------
    Index into array that is closest to val

    """
    return (np.abs(array-val)).argmin()


def time_slice(spikes, t_start, t_stop):
    """Slices into spike list between a start and stop time.

    Parameters
    ----------
    spikes : list of lists
        Where each inner list contains the spike times for an
        individual neuron. And len(spikes) is the total number
        of neurons.
    t_start : int
        If None, takes the slice from the beginning of the spike
        times.
    t_stop : int
        If None, takes the slice from the end of the spike times.

    Returns
    -------
    sliced_spikes : list of lists
        Where each inner list contains the spike times of interest
        for an individual neuron.

    Raises
    ------
    AssertionError
    When that len(spikes) != len(sliced_spikes) (eg. the number
    of neurons stays the same.
    """
    if t_start is None:
        t_start = -np.inf
    if t_stop is None:
        t_stop = np.inf

    sliced_spikes = []
    for neuron_spikes in spikes:
        indices = (neuron_spikes >= t_start) & (neuron_spikes <= t_stop)
        sliced_spikes.append(neuron_spikes[indices])

    assert(len(spikes) == len(sliced_spikes))

    return sliced_spikes


def idx_in_pos(position, index):
    """ Indexes into position data.

    Keeps x, y, time consistent.

    Parameters
    ----------
    position : dict
        With x, y, time keys
    index : int

    Returns
    -------
    pos : dict
        With x, y, time keys

    """
    pos = dict()
    pos['x'] = position['x'][index]
    pos['y'] = position['y'][index]
    pos['time'] = position['time'][index]
    return pos


def get_sort_idx(tuning_curves):
    """ Finds indices to sort neurons by max firing in tuning curve.

    Parameters
    ----------
    tuning_curves : list of lists
        Where each inner list is the tuning curves for an individual
        neuron.

    Returns
    -------
    sorted_idx : list
        List of integers that correspond to the neuron in sorted order.

    """
    tc_max_loc = []
    for i, neuron_tc in enumerate(tuning_curves):
        tc_max_loc.append((i, np.where(neuron_tc == np.max(neuron_tc))[0][0]))
    sorted_by_tc = sorted(tc_max_loc, key=lambda x: x[1])

    sorted_idx = []
    for idx in sorted_by_tc:
        sorted_idx.append(idx[0])

    return sorted_idx


# Adapted from mpl_toolkits.axes_grid2
class AnchoredScaleBar(AnchoredOffsetbox):
    def __init__(self, transform, sizex=0, sizey=0, labelx=None, labely=None,
                 loc=4, pad=0.1, borderpad=0.1, sep=2, prop=None, **kwargs):
        """
        Modified, draw a horizontal and/or vertical  bar with the size in data coordinate
        of the give axes. A label will be drawn underneath (center-aligned).

        Parameters
        ----------
        transform : the coordinate frame (typically axes.transData)
        sizex, sizey : width of x,y bar, in data units. 0 to omit
        labelx, labely : labels for x,y bars; None to omit
        loc : position in containing axes
        pad, borderpad : padding, in fraction of the legend font size (or prop)
        sep : separation between labels and bars in points.
        **kwargs : additional arguments passed to base class constructor
        """
        from matplotlib.lines import Line2D
        from matplotlib.pyplot import arrow
        from matplotlib.text import Text
        from matplotlib.offsetbox import AuxTransformBox
        bars = AuxTransformBox(transform)
        inv = transform.inverted()
        pixelxy = inv.transform((1, 1)) - inv.transform((0, 0))

        if sizex:
            barx = Line2D([sizex, 0], [0, 0], transform=transform, color='k')
            bars.add_artist(barx)

        if sizey:
            bary = Line2D([0, 0], [0, sizey], transform=transform, color='k')
            bars.add_artist(bary)

        if sizex and labelx:
            textx = Text(text=labelx, x=sizex/2.0, y=-5*pixelxy[1], ha='center', va='top')
            bars.add_artist(textx)

        if sizey and labely:
            texty = Text(text=labely, rotation='vertical', y=sizey/2.0, x=-2*pixelxy[0],
                         va='center', ha='right')
            bars.add_artist(texty)

        AnchoredOffsetbox.__init__(self, loc=loc, pad=pad, borderpad=borderpad,
                                       child=bars, prop=prop, frameon=False, **kwargs)

def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, **kwargs):
    """ Add scalebars to axes
    Adds a set of scale bars to *ax*, matching the size to the ticks of the
    plot and optionally hiding the x and y axes

    Parameters
    ----------
    ax : the axis to attach ticks to
    matchx, matchy : if True, set size of scale bars to spacing between ticks
                    if False, size should be set using sizex and sizey params
    hidex, hidey : if True, hide x-axis and y-axis of parent
    **kwargs : additional arguments passed to AnchoredScaleBars

    Returns created scalebar object
    """
    def find_loc(axis):
        loc = axis.get_majorticklocs()
        return len(loc)>1 and (loc[1] - loc[0])

    if matchx:
        kwargs['sizex'] = find_loc(ax.xaxis)
#         kwargs['labelx'] = str(kwargs['sizex'])
        kwargs['labelx'] = str(int(kwargs['sizex']*1000)) + ' ms'

    if matchy:
        kwargs['sizey'] = find_loc(ax.yaxis)
        kwargs['labely'] = str(kwargs['sizey'])

    scalebar = AnchoredScaleBar(ax.transData, **kwargs)
    ax.add_artist(scalebar)

    return scalebar
