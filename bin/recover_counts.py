#!/usr/bin/env python
"""
Reads in:
    ph*-eh*-accessibles: accessible states, energy, and counts

Writes out:
    ph*-eh*-accessibles.recovered: analytically recovered counts
"""
import os
import numpy as np
from pymcce import *

b = -KCAL2KT/(env.prm["MONTE_T"]/ROOMT)

def recover_counts(fn):
    lines = open(fn).readlines()
    n = len(lines)

    states = []
    Es = np.zeros(n)
    counts = np.zeros(n, dtype=int)
    for i in range(n):
        line = lines[i]
        state, value = line.split(":")
        E, count = [float(x) for x in value.split(",")]
        Es[i] = E
        counts[i] = count
        states.append(state)

    tared_Es = Es - Es[0]
    occ_state = np.exp(b*tared_Es)
    total_occ = np.sum(occ_state)
    occ_state_norm = occ_state / total_occ
    total_counts = np.sum(counts)
    recovered_counts = occ_state_norm * total_counts


    out_fn = "%s.recovered" % fn
    fh = open(out_fn, "w")
    for i in range(n):
        line = "%s:%.2f,%d,%d\n" % (states[i], Es[i], counts[i], recovered_counts[i])
        fh.write(line)
    fh.close()


    return

if __name__ == "__main__":
    folder = "microstates"

    # compose file names to read
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith("-accessibles")]


    for fn in files:
        print("Computing analytical recovered counts for %s" % fn)
        recover_counts(os.path.join(folder, fn))

