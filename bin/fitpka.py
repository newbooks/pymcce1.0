#!/usr/bin/env python
"""
Reads in:
    ph*-eh*-accessibles: accessible states, energy, and counts

Writes out:
    fort.38: occupancy table from sampling
    fort.38.recovery: occupancy table from analytical solution
    fort.38.xts: occupancy table with entropy correction
    sumcrg: Total charge table from sampling
    sumcrg.recovery: total charge table from analytical solution
    sumcrg.xts: total charge table with entropy correction
"""

from pymcce import *

if __name__ == "__main__":
    # compose file names to read
    folder = env.mc_states
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and f.endswith(".accessibles")]
    print(files)