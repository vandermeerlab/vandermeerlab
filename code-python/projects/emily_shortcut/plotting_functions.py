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


def raster_plot(spikes, colour='k'):
    location = 1
    for neuron in spikes:
        if len(neuron) > 0:
            plt.plot(neuron, np.ones(len(neuron))+location, '|', color=colour, ms=4, mew=1)
            location += 1
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron number')
    sns.despine()
    plt.ylim(0, location+1)


def plot_sorted_tc(sorted_tc, filepath):
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(sorted_tc, cmap='YlGn')
    plt.ylim(0, len(sorted_tc))
    plt.xlim(0, len(sorted_tc[0]))
    plt.ylabel('Neuron number')
    plt.xlabel('Location (cm)')
    sns.despine()
    #     plt.show()
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()


# Behavior
def plot_intersects(zone):
    for intersect in zone:
        plt.plot(intersect.exterior.xy[0], intersect.exterior.xy[1], 'b', lw=1)


def plot_zone(zone):
    plt.plot(zone.exterior.xy[0], zone.exterior.xy[1], 'b', lw=1)


def plot_bydurations(durations, filepath):
    ax = sns.boxplot(data=[durations['u'], durations['shortcut'], durations['novel']])
    sns.color_palette("hls", 18)
    ax.set(xticklabels=['U', 'Shortcut', 'Novel'])
    plt.ylabel('Duration of trial (s)')
    plt.xlabel('sessions=' + str(durations['num_sessions']))
    plt.ylim(0, 140)
    sns.despine()
    #     plt.show()
    plt.savefig(os.path.join(filepath, 'shortcut_behavior_durations.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_proportions(us, shortcuts, novels, filepath):
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

    plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join(filepath, 'shortcut_behaviour_proportions.png'), dpi=300, bbox_inches='tight')
    plt.close('all')


def plot_bytrial(togethers, filepath, min_length=30):
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
    # plt.show()
    plt.savefig(os.path.join(filepath, 'shortcut_behaviour_bytrial.png'), dpi=300, bbox_inches='tight')
    plt.close()


# Place fields
def plot_fields(heatmaps, pos, savepath, num, plot_log=True, num_bins=100, vmax=None):
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
    plt.text(2, 6, r'n=' + str(num), fontsize=15)
    # plt.show()
    plt.savefig(savepath, dpi=300)
    plt.close()


# Co-occurrence
def plot_cooccur(probs, filename):
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

    plt.tight_layout()
    # plt.show()
    print(filename)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


def plot_solo_coprob(probabilities, filepath, metric, title, ylabel):
    # Example function call:
    # plot_ind_coprob(probs, output_filepath, metric='p4', title='Co-activation above chance levels (p4)',
    # ylabel='SWR co-activation z-scored')
    ind = np.arange(3)
    width = 0.8
    colours = ['#5975a4', '#5f9e6e', '#b55d5f']
    fig, ax = plt.subplots()
    ax.bar(ind, [np.nanmean(probabilities['u'][metric]),
                 np.nanmean(probabilities['shortcut'][metric]),
                 np.nanmean(probabilities['novel'][metric])],
           width, color=colours)

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(ind + width * 0.5)
    ax.set_xticklabels(('U', 'Shortcut', 'Novel'))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.show()
    # plt.savefig(os.path.join(filepath, 'cooccur_', metric, '.png'), dpi=300, bbox_inches='tight')
    # plt.close()


def plot_swrs(csc, swr_idx, savepath, row=10, col=8, buffer=20):
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

        # plt.show()
        plt.savefig(savepath + str(fig + 1) + '.png', dpi=300)
        plt.close('all')
