from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import scipy.stats as stats
import seaborn as sns
import os

from behavior_functions import bytrial_counts, summary_bytrial

sns.set_style('white')
sns.set_style('ticks')


def raster_plot(spikes, savepath, savefig=False):
    """Plots raster plot of spikes from multiple neurons.

    Parameters
    ----------
    spikes : list of np.arrays
        Where each inner array contains the spike times (floats) for an individual
        neuron.
    savepath : str
        Location and filename for the saved plot.
    savefig : boolean
        Default is False show the plot without saving it. True and will save the
        plot to the specified location.

    """
    location = 1
    for neuron in spikes:
        if len(neuron) > 0:
            plt.plot(neuron, np.ones(len(neuron))+location, '|', color='k', ms=4, mew=1)
            location += 1
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron number')
    sns.despine()
    plt.ylim(0, location+1)

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_sorted_tc(sorted_tc, savepath, savefig=True):
    """Plots sorted tuning curves from a single trajectory for a session.

    Parameters
    ----------
    sorted_tc : list of lists
        Where each inner list contains the tuning curve (floats) for an individual
        neuron.
    savepath : str
        Location and filename for the saved plot.
    savefig : boolean
        Default is True and will save the plot to the specified location.
        False shows with plot without saving it.

    """
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(sorted_tc, cmap='YlGn')
    plt.ylim(0, len(sorted_tc))
    plt.xlim(0, len(sorted_tc[0]))
    plt.ylabel('Neuron number')
    plt.xlabel('Location (cm)')
    sns.despine()

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


# Behavior
def plot_intersects(zone):
    """Plots intersections of zones of ideal trajectory

    Parameters
    ----------
    zone : Shapely's Polygon object

    """
    for intersect in zone:
        plt.plot(intersect.exterior.xy[0], intersect.exterior.xy[1], 'b', lw=1)


def plot_zone(zone):
    """Plots zone of ideal trajectory

        Parameters
        ----------
        zone : Shapely's Polygon object

        """
    plt.plot(zone.exterior.xy[0], zone.exterior.xy[1], 'b', lw=1)


def plot_bydurations(durations, savepath, savefig=True):
    """Plots duration for each trial separated by trajectories. Behavior only.

        Parameters
        ----------
        durations : dict
            With u, shortcut, novel, num_sessions as keys.
            Each value is a list of durations (float) for a each session.
        savepath : str
            Location and filename for the saved plot.
        savefig : boolean
            Default is True and will save the plot to the specified location. False
            shows with plot without saving it.

        """
    ax = sns.boxplot(data=[durations['u'], durations['shortcut'], durations['novel']])
    sns.color_palette("hls", 18)
    ax.set(xticklabels=['U', 'Shortcut', 'Novel'])
    plt.ylabel('Duration of trial (s)')
    plt.xlabel('sessions=' + str(durations['num_sessions']))
    plt.ylim(0, 140)
    sns.despine()

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_proportions(us, shortcuts, novels, savepath, savefig=True):
    """Plots proportion of each trajectory taken. Behavior only.

        Parameters
        ----------
        us : list of floats
            Proportion along the u trajectory for each session.
            len(us) == num_sessions evaluated
        shortcuts : list of floats
            Proportion along the shortcut trajectory for each session.
            len(shortcut) == num_sessions evaluated
        novels : list of floats
            Proportion along the novel trajectory for each session.
            len(novel) == num_sessions evaluated
        savepath : str
            Location and filename for the saved plot.
        savefig : boolean
            Default is True and will save the plot to the specified location. False
            shows with plot without saving it.

        """
    all_us = np.mean(us)
    us_sem = stats.sem(us)
    all_shortcuts = np.mean(shortcuts)
    shortcuts_sem = stats.sem(shortcuts)
    all_novels = np.mean(novels)
    novels_sem = stats.sem(novels)

    n_groups = list(range(3))

    colour = ['#5975a4', '#5f9e6e', '#b55d5f']

    data = [all_us, all_shortcuts, all_novels]
    sems = [us_sem, shortcuts_sem, novels_sem]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in list(range(len(data))):
        ax.bar(n_groups[i], data[i], align='center',
               yerr=sems[i], color=colour[i], ecolor='#525252')

    plt.xlabel('(sessions=' + str(len(us)) + ')')
    plt.ylabel('Proportion of trials')
    sns.despine()
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.xticks(n_groups, ['U', 'Shortcut', 'Novel'])

    # plt.tight_layout()
    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_bytrial(togethers, savepath, min_length=30, savefig=True):
    """Plots choice of trajectory by trial. Behavior only.

        Parameters
        ----------
        togethers : list
        savepath : str
            Location and filename for the saved plot.
        min_length = int
            This is the number of trials to be considered.
            The default is set to 30 (Eg. trials 1-30 are considered).
        savefig : boolean
            Default is True and will save the plot to the specified location.
            False shows with plot without saving it.

        """
    bytrial = bytrial_counts(togethers, min_length)

    means, sems = summary_bytrial(bytrial, min_length)

    trials = list(range(min_length))

    colours = dict(u='#5975a4', shortcut='#5f9e6e', novel='#b55d5f')
    labels = dict(u='Full U', shortcut='Shortcut', novel='Novel')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for path in means:
        ax.plot(trials, means[path], color=colours[path], label=labels[path], marker='o', lw=2)
        ax.fill_between(trials, np.array(means[path]) - np.array(sems[path]),
                        np.array(means[path]) + np.array(sems[path]),
                        color=colours[path], interpolate=True, alpha=0.3)
    plt.ylabel('Proportion of trials')
    plt.xlabel('Trial number (sessions=' + str(len(togethers)) + ')')
    sns.despine()
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.legend(loc=1, prop={'size': 10})

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


# Place fields
def plot_fields(heatmaps, pos, num_neurons, savepath, savefig=True, plot_log=True, num_bins=100, vmax=None):
    """Plots spiking in 2D position space.

        Parameters
        ----------
        heatmaps : dict of lists
            Where the key is the neuron number and the value is the heatmap for
            that individual neuron.
        pos : dict
            With time(float), x(float), y(float) as keys
        savepath : str
            Location and filename for the saved plot.
        num_neurons = int
            Number of neurons used for this plot.
        savefig : boolean
            Default is True and will save the plot to the specified location.
            False shows with plot without saving it.
        plot_log : boolean
            Range for plot in log scale. Allows for better visualization of
            areas with large amount of spikes. Default is set to True.
        num_bins : int
            Number of bins the 2D space is divided into.
            Default is set to 100.
        vmax : float, optional
            Default is set to None. Used to scale the heatmaps for different
            neurons or sessions to better directly compare.

        """
    plt.figure()
    plt.plot(pos['x'], pos['y'], 'k.', ms=0.2)
    xedges = np.linspace(np.min(pos['x'])-2, np.max(pos['x'])+2, num_bins+1)
    yedges = np.linspace(np.min(pos['y'])-2, np.max(pos['y'])+2, num_bins+1)
    xx, yy = np.meshgrid(xedges, yedges)

    if plot_log:
        pp = plt.pcolormesh(xx, yy, heatmaps, norm=SymLogNorm(linthresh=1.0, vmax=vmax), cmap='YlGn')
    else:
        pp = plt.pcolormesh(xx, yy, heatmaps, vmax=vmax, cmap='YlGn')
    print(np.max(heatmaps))
    plt.colorbar(pp)
    plt.axis('off')
    plt.text(2, 6, r'n=' + str(num_neurons), fontsize=15)

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


# Co-occurrence
def plot_cooccur(probs, savepath, savefig=True):
    """Plots co-occurrence probabilities from p0, p2, p3, p4.

        Parameters
        ----------
        probs : dict
            Where the keys are p0, p1, p2, p3, p4, p5. P1 through p4 are
            np.arrays of length num_neurons
        savepath : str
            Location and filename for the saved plot.
        savefig : boolean
            Default is True and will save the plot to the specified location.
            False shows with plot without saving it.

        Notes
        -----
        p0 : probability (fraction of time bins) each neuron is active.
        p2 : Observed conditional probability
            .. math:: p(x|y)
        p3 : Observed co-occurrence (joint) probability
            .. math:: p(x,y)
        p4 : z-score of p3 against shuffled data

        """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ind = np.arange(3)
    width = 0.8
    colours = ['#5975a4', '#5f9e6e', '#b55d5f']

    ax1.bar(ind, [np.nanmean(probs['u']['p0']), np.nanmean(probs['shortcut']['p0']), np.nanmean(probs['novel']['p0'])],
            width, color=colours)
    ax1.set_ylabel('Proportion of SWRs active')
    ax1.set_title('Probability that a given neuron is active (p0)')
    ax1.set_xticks(ind + width*0.5)
    ax1.set_xticklabels(('U', 'Shortcut', 'Novel'))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    ax2.bar(ind, [np.nanmean(probs['u']['p2']), np.nanmean(probs['shortcut']['p2']), np.nanmean(probs['novel']['p2'])],
            width, color=colours)
    ax2.set_ylabel('??? Some sort of probability')
    ax2.set_title('Observed conditional probabilities (p2)')
    ax2.set_xticks(ind + width*0.5)
    ax2.set_xticklabels(('U', 'Shortcut', 'Novel'))
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()

    ax3.bar(ind, [np.nanmean(probs['u']['p3']), np.nanmean(probs['shortcut']['p3']), np.nanmean(probs['novel']['p3'])],
            width, color=colours)
    ax3.set_ylabel('Cell pair joint probability')
    ax3.set_title('Observed co-activity (p3)')
    ax3.set_xticks(ind + width*0.5)
    ax3.set_xticklabels(('U', 'Shortcut', 'Novel'))
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()

    ax4.bar(ind, [np.nanmean(probs['u']['p4']), np.nanmean(probs['shortcut']['p4']), np.nanmean(probs['novel']['p4'])],
            width, color=colours)
    ax4.set_ylabel('SWR co-activation z-scored')
    ax4.set_title('Co-activation above chance levels (p4)')
    ax4.set_xticks(ind + width*0.5)
    ax4.set_xticklabels(('U', 'Shortcut', 'Novel'))
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.get_xaxis().tick_bottom()
    ax4.get_yaxis().tick_left()

    # plt.tight_layout()

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_solo_coprob(probs, savepath, metric, title, ylabel, savefig=False):
    """Plots co-occurrence probabilities from p0, p2, p3, p4.

        Parameters
        ----------
        probs : dict
            Where the keys are p0, p1, p2, p3, p4, p5. P1 through p4 are
            np.arrays of length num_neurons
        savepath : str
            Location and filename for the saved plot.
        metric : str
            One of the keys in probs.
        title : str
        ylabel : str
        savefig : boolean
            Default is False and will show the plot without saving it.
            True will save the plot to the specified location.

        Notes
        -----
        p0 : probability (fraction of time bins) each neuron is active.
        p1 : expected co-occurrence under independence assumption
            .. math:: p(x,y) = p(x) * p(y)
        p2 : Observed conditional probability
            .. math:: p(x|y)
        p3 : Observed co-occurrence (joint) probability
            .. math:: p(x,y)
        p4 : z-score of p3 against shuffled data
        p5 : Observed co-occurrence (joint) probability of shuffled data
            .. math:: p(x,y)

        """
    ind = np.arange(3)
    width = 0.8
    colours = ['#5975a4', '#5f9e6e', '#b55d5f']
    fig, ax = plt.subplots()
    ax.bar(ind, [np.nanmean(probs['u'][metric]),
                 np.nanmean(probs['shortcut'][metric]),
                 np.nanmean(probs['novel'][metric])],
           width, color=colours)

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(ind + width * 0.5)
    ax.set_xticklabels(('U', 'Shortcut', 'Novel'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    if savefig:
        plt.savefig(savepath, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_swrs(csc, swr_idx, saveloc, row=10, col=8, buffer=20, savefig=True):
    """Plots all local field potentials (LFP) around sharp-wave ripple (SWR) times
        for each given SWR.

        Parameters
        ----------
        csc : dict
            With time(np.array), data(np.array) as keys
        swr_idx : dict
            With start(int), stop(int) as keys
        saveloc : str
            Location and filename for the saved plot. Do not add '.png', it is
            added here to include the multiple figures into the filename.
        row = int
            Number of rows in each figure.
        col : int
            Number of columns in each figure.
        buffer : int
            Amount of LFP shown on either side of SWR.
        savefig : boolean
            Default is True and will save the plot to the specified location.
            False shows with plot without saving it.

        """
    plots_per_fig = row * col
    num_figures = range(int(np.ceil(len(swr_idx['start']) / plots_per_fig)))

    for fig in num_figures:
        print('figure', fig, 'of', np.max(list(num_figures)))
        plt.figure(fig)

        stop_idx = plots_per_fig * (fig + 1)
        start_idx = stop_idx - plots_per_fig
        if stop_idx > len(swr_idx['start']):
            stop_idx = len(swr_idx['start']) + 1

        for i, (start, stop) in enumerate(
                zip(swr_idx['start'][start_idx:stop_idx], swr_idx['stop'][start_idx:stop_idx])):
            plt.subplot(row, col, i + 1)

            plt.plot(csc['time'][start - buffer:stop + buffer], csc['data'][start - buffer:stop + buffer], 'k')
            plt.plot(csc['time'][start:stop], csc['data'][start:stop], 'r')

            plt.axis('off')

        if savefig:
            plt.savefig(saveloc + str(fig + 1) + '.png', dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
