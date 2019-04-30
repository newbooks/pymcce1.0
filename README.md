# pymcce1.0
Python MCCE

Output:
fort.38 native, recovered and TS corrected
titration curve
accessible states, steps, total states
sampling energy trace : traceenergy.py
states vs energy
occupancy vs energy

valid_state(state):
    each conf in state is in free_residues
    each res in free_residues has one and only one conf in state
    This makes sure the state is free residues only and garauntees correct energy

To calculate system total energy:
    fixed conformers impact other conformers as mfe self energy
    sum up self energy of each conformer with occupancy
    subtract pairwise interaction between fixed conformers, because they were conunted twice in previous steps.
