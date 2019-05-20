#!/usr/bin/env python
import random
import math
import sys

N = int(sys.argv[1])

N_eq = 1000
N_state = 100
state_Es = [i*0.1 for i in range(N)]
state = random.randrange(N_state)
counters = [0] * N_state

def mc(state):
    new = random.randrange(N_state)
    dE = state_Es[new] - state_Es[state]
    if dE < 0:
        i_state = new
    elif random.random() < math.exp(-dE):
        i_state = new
    else:
        i_state = state

    return i_state

# equlibration
for i in range(N_eq):
    state = mc(state)

# sampling
for i in range(N):
    state = mc(state)
    counters[state] += 1

print("E    Count")
for i in range(N_state):
    i_state = N_state - i -1
    print("%4.1f %4d" % (state_Es[i_state], counters[i_state]))