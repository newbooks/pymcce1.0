#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import style
import pandas as pd
import os
import sys


def show_stateP(fname):
    n_intervals = 50

    style.use('ggplot')
    data = pd.read_csv(fname, delimiter = ':')
    value = data.iloc[:,1]
    c = value.str.split(",", expand=True)
    E = c.iloc[:,0].values.astype(np.float)
    MC_counts = c.iloc[:,1].values.astype(np.int)
    Re_counts = c.iloc[:,2].values.astype(np.int)


    fig, windows = plt.subplots(2, sharex = True)
    windows[0].plot(E, MC_counts, label="Monte Carlo")
    windows[0].legend(loc="upper right")
    windows[0].plot(E, Re_counts, label="Analytical")
    windows[0].legend(loc="upper right")
    windows[0].set_title("counts of individual state")

    Emin = E[0]
    interval = (E[-1]-E[0])/n_intervals    # 100 intervals along energy
    x = [Emin + 0.5 * interval + interval * i for i in range(n_intervals)]   # mid point of intervals as x ticks
    y1 = np.zeros(n_intervals)
    y2 = np.zeros(n_intervals)
    for i in range(len(E)):
        j = int((E[i] - Emin)/interval)
        if j >= n_intervals:
            j = n_intervals - 1   # last may sit on the last tick causing index overflow
        y1[j] += MC_counts[i]
        y2[j] += Re_counts[i]


    bar_width = interval/2
    #windows[1].plot(x, y1, label="Monte Carlo")
    windows[1].bar(x, y1, bar_width, label="Monte Carlo")
    windows[1].legend(loc="upper right")
    #windows[1].plot(x, y2, label="Analytical")
    windows[1].bar(x+bar_width, y2, bar_width, label="Analytical")
    windows[1].legend(loc="upper right")
    windows[1].set_title("counts of energy band")

    plt.show()

    return


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Specify a file of accessible states.")
        print("Example: visual_pvse.py microstates/ph0.0-eh0-accessibles.recovered")
        sys.exit()

    fname = sys.argv[1]
    show_stateP(fname)