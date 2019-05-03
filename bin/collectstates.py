#!/usr/bin/env python
"""
Collect unique states from microstates directory.
It reads in:
    all ms.gz files
It writes out:
    ph*-eh*-converge: divide the states into 6 runs x 20 groups, show average energy abd stdev of each
    ph*-eh*-states:   after discarding a percentage of eq runss, collect accessible states, energy, and counts
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

class State_stat:
    def __init__(self, E):
        self.E = E
        self.counter = 1
        return


def traceenergy(fname):
    with gzip.open(fname, "rb") as fh:
        lines = fh.readlines()

    fname_base = os.path.splitext(fname)[0]
    fname_energy = fname_base + ".energy"

    # head line, T, ph and eh
    line = lines.pop(0).decode()
    fields = line.strip().split(",")
    mc_parm = {}
    for field in fields:
        key, value = field.split("=")
        mc_parm[key.strip()] = float(value)

    T = mc_parm["T"]
    ph = mc_parm["ph"]
    eh = mc_parm["eh"]

    prot.update_energy(T=T, ph=ph, eh=eh)
    line = lines.pop(0).decode()
    _, state_str = line.split(":")
    if state_str:
        state = [int(ic) for ic in state_str.split()]
    else:
        print("No initial state found. Quitting ...")
        return

    if validate_state(prot, state):
        E = get_state_energy(prot, state)
    else:
        print("Quitting ...")
        return

    fh = open(fname_energy, "w")

    state = set(state)
    for line in lines:
        line = line.strip().decode()
        if line:
            delta = line.split(",")
            for ic in delta:
                ic = int(ic)
                if ic < 0:
                    ic = -ic
                    state = state - {ic}
                else:
                    state.add(ic)

            E = get_state_energy(prot, list(state))

        fh.write("%.3f\n" % E)

    fh.close()




    return



def collect_one(c, t):
    print("collecting microstates at %s and throw_away = %.2f%%. " % (c, t*100))
    # get files at this condition
    folder = "microstates"
    pattern = os.path.join(folder, "%s-run*.ms.gz" % c)
    files = glob.glob(pattern)
    files.sort()

    all_states = {}

    for f in files:
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
                all_states[state_tup] = State_stat(E)
            else:
                E = get_state_energy(prot, state_arr)
                all_states[state_tup] = State_stat(E)

        print(len(lines[n_skip:]))
        print(len(all_states))

    return



def collect(throwaway):
    # compose file names to read
    folder = "microstates"
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
