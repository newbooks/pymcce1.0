#!/usr/bin/env python
"""
Collect unique states from microstates directory.
It reads in:
    all ms.gz files
It writes out:
    ph*-eh*-converge: divide the states into 6 runs x 20 groups, show average energy abd stdev of each
    ph*-eh*-accessibles:   after discarding a percentage of eq runss, collect accessible states, energy,
    and counts
    fort.38: occupancy table from sampling
    fort.38.recovery: occupancy table from analytical solution
    fort.38.xts: occupancy table with entropy correction
    sumcrg: Total charge table from sampling
    sumcrg.recovery: total charge table from analytical solution
    sumcrg.xts: total charge table with entropy correction
"""

import sys
from pymcce import *
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
    # get files at this condition
    folder = env.mc_states
    pattern = os.path.join(folder, "%s-run*.ms.gz" % c)
    files = glob.glob(pattern)
    files.sort()

    all_states = {}

    std_stat = []
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

        prot.update_energy(T=T, ph=ph, eh=eh)
        line = lines.pop(0)
        _, state_str = line.split(":")
        if state_str:
            state = [int(ic) for ic in state_str.split()]
        else:
            print("   No initial state found. Quitting ...")
            continue

        if not validate_state(prot, state):
            print("   Quiting ...")
            continue

        # Now we have initial state in state[], and the rest state deltas in lines[]
        # skip t lines
        n_lines = len(lines)
        n_skip = int(t * n_lines)
        state = set(state)
        for line in lines[:n_skip]:
            line = line.strip()
            if line:
                delta = line.split(",")
                for ic in delta:
                    ic = int(ic)
                    if ic < 0:
                        ic = -ic
                        state = state - {ic}
                    else:
                        state.add(ic)

        # collect the rest states
        state_arr = list(state)
        state_arr.sort()
        state_tup = tuple(state_arr)
        E = get_state_energy(prot, state_arr)   # initial state

        n_record = len(lines[n_skip:])
        n_segment = int(n_record/20)

        counter_line = 0
        Es = np.zeros(n_segment)
        for line in lines[n_skip:]:
            line = line.strip()
            if line:
                delta = line.split(",")
                for ic in delta:
                    ic = int(ic)
                    if ic < 0:
                        ic = -ic
                        state = state - {ic}
                    else:
                        state.add(ic)
                # got an update
                state_arr = list(state)
                state_arr.sort()
                state_tup = tuple(state_arr)

            # save it to database
            if state_tup in all_states:
                all_states[state_tup].counter += 1
            else:
                E = get_state_energy(prot, state_arr)
                all_states[state_tup] = State_stat(E)

            # stdev
            counter_line += 1
            Es[counter_line-1] = E
            if counter_line >= n_segment:
                std_stat_f.append((Es.mean(), Es.std()))
                counter_line = 0
        std_stat.append(std_stat_f)


    fn_stats = "%s/%s-accessibles_stats" % (folder, c)
    out_lines = ["Segments  %12s\n" % (" ".join([f.split("-")[-1].split(".")[0] for f in files]))]
    for i in range(20):
        stat_str = ["%5.2f|%5.2f" % (std_stat[j][i][0], std_stat[j][i][1]) for j in range(len(files))]
        out_lines.append("%-12d %s\n" % (i+1, " ".join(stat_str)))
    open(fn_stats, "w").writelines(out_lines)

    fn_accessibles = "%s/%s-accessibles" % (folder, c)
    accessibles = sorted(all_states.items(), key=lambda kv:kv[1].E)
    out_lines = []
    for rs in accessibles:
        out_lines.append("%s:%.2f, %d\n" %(rs[0], rs[1].E, rs[1].counter))
    open(fn_accessibles, "w").writelines(out_lines)

    return



def collect(throwaway):
    # compose file names to read
    folder = env.mc_states
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith(".ms.gz")]
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
    prot = MC_Protein()
    if len(sys.argv) < 2:
        print("collectstates.py cutoff_percentage (default 0.1)")
        print("          cutoff_percentage is the fraction of initial microstates to throw away")
        sys.exit()
    throwaway = float(sys.argv[1])
    collect(throwaway)
