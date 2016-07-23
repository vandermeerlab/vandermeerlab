import numpy as np
import scipy.stats as stats


def bytrial_counts(togethers, min_length):
    bytrial = dict(u=[], shortcut=[], novel=[])
    for trial in list(range(min_length)):
        for key in bytrial:
            bytrial[key].append([])

    for session in togethers:
        for trial in list(range(min_length)):
            for key in bytrial:
                if session[trial][1] == key:
                    bytrial[key][trial].append(1)
                else:
                    bytrial[key][trial].append(0)
    return bytrial


def summary_bytrial(bytrial, min_length):
    means = dict(u=[], shortcut=[], novel=[])
    sems = dict(u=[], shortcut=[], novel=[])
    for trial in list(range(min_length)):
        for key in bytrial:
            means[key].append(np.mean(bytrial[key][trial]))
            sems[key].append(stats.sem(bytrial[key][trial]))
    return means, sems
