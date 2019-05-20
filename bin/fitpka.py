#!/usr/bin/env python
"""
Reads in:
    ph*-eh*-accessibles: accessible states, energy, and counts

Writes out:
    ph*-eh*-accessibles.recovered: analytically recovered counts
    fort.38: occupancy table from sampling
    fort.38.recovery: occupancy table from analytical solution
    fort.38.xts: occupancy table with entropy correction
    sumcrg: Total charge table from sampling
    sumcrg.recovery: total charge table from analytical solution
    sumcrg.xts: total charge table with entropy correction
"""
import os
import sys

PH2KCAL = 1.364

def titration_range(files):
    ph_range = []
    eh_range = []

    for fn in files:
        fields = fn.split("-")
        ph = fields[0].split("ph")[1]
        eh = fields[1].split("eh")[1]
        if ph not in ph_range:
            ph_range.append(ph)
        if eh not in eh_range:
            eh_range.append(eh)

    ph_range = [float(x) for x in ph_range]
    eh_range = [float(x) for x in eh_range]
    ph_range.sort()
    eh_range.sort()
    return ph_range, eh_range

def recover_counts(fn):
    lines = open(fn).readlines()

    return

if __name__ == "__main__":
    # compose file names to read
    folder = "microstates"
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith("-accessibles")]

    ph_range, eh_range = titration_range(files)
    if len(ph_range) > 1 and len(eh_range) == 1:
        titration_type = "ph"
        points = ph_range
    elif len(ph_range) == 1 and len(eh_range) > 1:
        titration_type = "eh"
        points = eh_range
    elif len(ph_range) == 1 and len(eh_range) == 1:
        titration_type = "single"
    else:
        titration_type = "multi"

    for fn in files:
        print("Computing analytical recovered counts)")
        recover_counts(fn)

    ph_range_str = ", ".join(["%.1f" %x for x in ph_range])
    eh_range_str = ", ".join(["%.0f" %x for x in eh_range])
    print("Titration type is \"%s\" at pH [%s] and eh [%s]" % (titration_type, ph_range_str, eh_range_str))
    if titration_type == "multi":
        print("Skipping curve fitting")
        sys.exit()

