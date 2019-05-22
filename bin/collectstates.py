#!/usr/bin/env python
"""
Collect unique states from microstates directory.
It reads in:
    all ms.gz files
It writes out:
    ph*-eh*-accessibles.stats: divide the states into 6 runs x 20 groups, show average energy abd stdev of each
    ph*-eh*-accessibles:  after discarding a percentage of eq runs, collect accessible states, energy,
    and counts
"""

import sys
import gzip
import os
import glob
import numpy as np

class State_stat:
    def __init__(self, E):
        self.E = E
        self.counter = 1
        return


def collect_one(c, t):
    print("collecting microstates at %s and throw_away = %.2f%%. " % (c, t*100))
    visited = 0
    # get files at this condition
    folder = "microstates"
    pattern = os.path.join(folder, "%s-run*.ms.gz" % c)
    files = glob.glob(pattern)
    files.sort()

    all_states = {}

    std_stat = []
    number_of_acc = []
    for f in files:
        std_stat_f = []
        print("   Processing file %s" % f)
        with gzip.open(f, "rb") as fh:
            lines = fh.readlines()
        lines = [x.decode() for x in lines]

        line = lines.pop(0)
        fields = line.strip().split(",")
        mc_parm = {}
        for field in fields:
            key, value = field.split("=")
            mc_parm[key.strip()] = float(value)

        T = mc_parm["T"]
        ph = mc_parm["ph"]
        eh = mc_parm["eh"]

        line = lines.pop(0)
        E_str, state_str = line.split(":")
        if state_str:
            state = [int(ic) for ic in state_str.split(",")]
            E = float(E_str)
        else:
            print("   No initial state found. Quitting ...")
            continue

        # Now we have initial state in state[], and the rest state deltas in lines[]
        # skip t lines
        n_lines = len(lines)
        n_skip = int(t * n_lines)
        state = set(state)
        for line in lines[:n_skip]:
            line = line.strip()
            if line:
                fields = line.split(":")
                E = float(fields[0])
                off_confs, onconfs = conf_delta(fields[1])
                state = state - off_confs
                state = state | onconfs

        # collect the rest states
        state_arr = list(state)
        state_arr.sort()
        state_tup = tuple(state_arr)

        n_record = len(lines[n_skip:])
        n_segment = int(n_record/20)

        counter_line = 0
        Es = np.zeros(n_segment)
        for line in lines[n_skip:]:
            line = line.strip()
            if line:
                fields = line.split(":")
                E = float(fields[0])
                off_confs, onconfs = conf_delta(fields[1])
                state = state - off_confs
                state = state | onconfs

                # got an update
                state_arr = list(state)
                state_arr.sort()
                state_tup = tuple(state_arr)

            # save it to database
            if state_tup in all_states:
                all_states[state_tup].counter += 1
            else:
                all_states[state_tup] = State_stat(E)

            # stdev
            counter_line += 1
            Es[counter_line-1] = E
            if counter_line >= n_segment:
                std_stat_f.append((Es.mean(), Es.std()))
                counter_line = 0

        visited += len(lines) - n_skip

        std_stat.append(std_stat_f)
        number_of_acc.append((len(all_states), visited))


    fn_stats = "%s/%s-accessibles.stats" % (folder, c)
    out_lines = ["Segments          %s\n" % ("               ".join([f.split("-")[-1].split(".")[0] for f in files]))]
    for i in range(20):
        stat_str = ["%9.2f/%-9.2f" % (std_stat[j][i][0], std_stat[j][i][1]) for j in range(len(files))]
        out_lines.append("%-12d %s\n" % (i+1, " ".join(stat_str)))

    stat_str = " ".join(["%8d/%-10d" % (n[0], n[1]) for n in number_of_acc])
    out_lines.append("#:            %s\n" % stat_str)

    open(fn_stats, "w").writelines(out_lines)

    fn_accessibles = "%s/%s-accessibles" % (folder, c)
    accessibles = sorted(all_states.items(), key=lambda kv:kv[1].E)
    out_lines = []
    for rs in accessibles:
        out_lines.append("%s:%.3f, %d\n" %(rs[0], rs[1].E, rs[1].counter))
    open(fn_accessibles, "w").writelines(out_lines)

    return


def conf_delta(line):
    off_confs = set()
    on_confs = set()
    delta = line.split(",")
    for ic_s in delta:
        if ic_s[0] == "-":
            ic = -int(ic_s)
            off_confs.add(ic)
        else:
            ic = int(ic_s)
            on_confs.add(ic)
    return off_confs, on_confs

def collect(throwaway):
    # compose file names to read
    folder = "microstates"
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith(".ms.gz")]
    files.sort()
    # analyze how many ph-eh
    titr_conditions = []
    for f in files:
        fields = f.split("-")
        c = "-".join([fields[0], fields[1]])
        if c not in titr_conditions:
            titr_conditions.append(c)

    print("")
    for c in titr_conditions:
        collect_one(c, throwaway)

    return


if __name__ == "__main__":
    if len(sys.argv) >1:
        throwaway = float(sys.argv[1])
    else:
        throwaway = 0.1
    collect(throwaway)
