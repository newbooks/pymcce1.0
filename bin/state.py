#!/usr/bin/env python

from pymcce import *

def analyze_state(fname):
    lines = open(fname).readlines()

    # head line, T, ph and eh
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
    print("Environment: pH = %.2f eh = %.0f Temperature = %.2f K" % (ph, eh, T))
    for line in lines:
        state = [int(ic) for ic in line.split(",")]

        if validate_state(prot, state):
            E = get_state_energy(prot, state)
            print("E = %.2f" % E)
        else:
            print("Not a valid state, Quitting ...")
            return

    return


if __name__ == "__main__":
    prot = MC_Protein()
    fname = "states"
    analyze_state(fname)
