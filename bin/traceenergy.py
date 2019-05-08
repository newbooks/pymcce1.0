#!/usr/bin/env python

import sys
from pymcce import *
import gzip
import os

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
            off_confs, onconfs = conf_delta(line)
            state = state - off_confs
            state = state | onconfs
            E = get_state_energy(prot, list(state))

        fh.write("%.3f\n" % E)

    fh.close()




    return


if __name__ == "__main__":
    prot = MC_Protein()
    fname = sys.argv[1]
    traceenergy(fname)
