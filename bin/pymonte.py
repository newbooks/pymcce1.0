#!/usr/bin/env python

print("Monte Carlo sampling:")

from pymcce import *
import time

if __name__ == "__main__":
    timerA = time.time()
    env.print_scaling()
    prot = MC_Protein()


    monte_t = env.prm["MONTE_T"]
    titration_type = env.prm["TITR_TYPE"].upper()
    init_ph = env.prm["TITR_PH0"]
    step_ph = env.prm["TITR_PHD"]
    init_eh = env.prm["TITR_EH0"]
    step_eh = env.prm["TITR_EHD"]
    steps = env.prm["TITR_STEPS"]

    total_states = 1
    for res in prot.free_residues:
        total_states *= len(res)

    timerB = time.time()
    print("   Done setting up MC in %1d seconds.\n" % (timerB - timerA))

    if total_states > env.prm["NSTATE_MAX"]:
        print("   Total states %d > threshold %d" % (total_states, env.prm["NSTATE_MAX"]))
    else:
        print("   Total states %d <= threshold %d" % (total_states, env.prm["NSTATE_MAX"]))
        print("   Please use analytical method ______ to analyze protein equilibrium0.")
        sys.exit()

    mc_prepdir()
    os.chdir(env.mc_states)

    for i in range(steps):
        # Set up pH and eh environment
        if titration_type == "PH":
            ph = init_ph + i * step_ph
            eh = init_eh
        elif titration_type == "EH":
            ph = init_ph
            eh = init_eh + i * step_eh
        else:
            print(
                "   Error: Titration type is %s. It has to be ph or eh in line (TITR_TYPE) in run.prm" % titration_type)
            sys.exit()


        mc_sample(prot, T=monte_t, ph=ph, eh=eh)

    os.chdir("../")