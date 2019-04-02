#!/usr/bin/env python

import sys
import os
import numpy as np

class Env:
    def __init__(self):
        # Hard coded values
        self.runprm= "run.prm"
        self.version = "PyMCCE 1.0"
        self.fn_conflist1 = "head1.lst"
        self.fn_conflist2 = "head2.lst"
        self.fn_conflist3 = "head3.lst"
        self.energy_table = "energies"
        self.prm = self.load_runprm()
        self.tpl = {}
        self.read_extra()
        return

    def load_runprm(self):
        float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD", "CLASH_DISTANCE",
                        "BIG_PAIRWISE", "MONTE_T", "MONTE_REDUCE"]
        int_values = ["TITR_STEPS", "MONTE_RUNS", "MONTE_TRACE", "MONTE_NITER", "MONTE_NEQ",
                      "MONTE_NSTART", "MONTE_FLIPS", "NSTATE_MAX"]
        prm = {}
        print("   Loading %s" % self.runprm)
        lines = open(self.runprm).readlines()
        # Sample line: "t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
        for line in lines:
            line = line.strip()
            line = line.split("#")[0]  # This cuts off everything after #
            left_p = line.rfind("(")
            right_p = line.rfind(")")
            if left_p > 0 and right_p > left_p + 1:
                key = line[left_p + 1:right_p]
                fields = line[:left_p].split()
                if len(fields) >= 1:
                    value = fields[0]
                    if key in float_values:
                        prm[key] = float(value)
                    elif key in int_values:
                        prm[key] = int(value)
                    else:
                        prm[key] = value
        return prm

    def print_runprm(self):
        for key in self.prm.keys():
            print("%-25s:%s" % (key, str(self.prm[key])))
        return

    def load_ftpl(self, file):
        """Load a tpl file."""
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        print("   Loading ftpl file %s" % file)
        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys = tuple(keys)

            value_string = fields[1].strip()
            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return


    def load_tpl(self, file):
        """Load a tpl file."""
        print("   Loading tpl file %s" % file)
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            if len(line) < 21:
                continue
            keys = [line[:9], line[9:14], line[15:19]]
            value_string = line[20:].strip()

            keys = [x for x in keys if x]
            keys = tuple(keys)

            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return


    def read_extra(self):
        """Read extra.tpl."""
        fname = self.prm["EXTRA"]

        print("   Extra tpl parameters in file %s" % fname)
        if os.path.isfile(fname):
            if fname[-5:] == ".ftpl":
                self.load_ftpl(fname)
            elif fname[-4:] == ".tpl ":
                self.load_tpl(fname)

        default_values_keys = [("SCALING", "VDW0"),
                               ("SCALING", "VDW1"),
                               ("SCALING", "VDW"),
                               ("SCALING", "TORS"),
                               ("SCALING", "ELE"),
                               ("SCALING", "DSOLV")]
        for element in default_values_keys:
            if element not in self.tpl:
                print("      Set to default: %s = 1.0" % ",".join(element))
                self.tpl[element] = 1.0

        return

    def print_scaling(self):
        """Print scaling factors."""
        # print self.param
        print("   Scaling factors:")
        print("   VDW0  = %.3f" % self.tpl[("SCALING", "VDW0")])
        print("   VDW1  = %.3f" % self.tpl[("SCALING", "VDW1")])
        print("   VDW   = %.3f" % self.tpl[("SCALING", "VDW")])
        print("   TORS  = %.3f" % self.tpl[("SCALING", "TORS")])
        print("   ELE   = %.3f" % self.tpl[("SCALING", "ELE")])
        print("   DSOLV = %.3f" % self.tpl[("SCALING", "DSOLV")])
        print("   Done\n")
        return


class Conformer:
    def __init__(self, fields):
        # from head3.lst
        self.iConf = int(fields[0])
        self.confname = fields[1]
        self.flag = fields[2].lower()
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = float(fields[5])
        self.pk0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9]) * env.param[("SCALING", "VDW0")]
        self.vdw1 = float(fields[10]) * env.param[("SCALING", "VDW1")]
        self.tors = float(fields[11]) * env.param[("SCALING", "TORS")]
        self.epol = float(fields[12]) * env.param[("SCALING", "ELE")]
        self.dsolv = float(fields[13]) * env.param[("SCALING", "DSOLV")]
        self.extra = float(fields[14])
        self.history = fields[15]
        # from MC entropy sampling
        self.entropy = 0.0  # -TS, will be calculated at entropy sampling
        self.acc_entropy = []
        # needed by MC process
        self.E_self = 0.0  # self energy in head3.lst
        self.E_self_mfe = 0.0  # self energy including pairwise contribution from fixed residues
        self.counter = 0  # MC counters
        self.mc_occ = 0.0  # MC calculated occ
        self.acc_occ = []  # a list of past mc_occ until rest, history is needed for test convergence
        self.occ_at_points = []  # occ at titration points
        return

class MC_Protein:
    """Monte Carlo Protein data structure."""
    def __init__(self):
        print("\n   Initializing Monte Carlo frame strcuture.")
        self.head3list = []
        self.confnames = []
        return

    def read_head3list(self):
        fname = env.fn_conflist3
        print("      Loading confomer self energy from %s" % fname)

        lines = open(fname).readlines()
        lines.pop(0)
        for line in lines:
            fields = line.split()
            if len(fields) >= 16:
                conf = Conformer(fields)
                self.head3list.append(conf)

        # validate
        self.confnames = [x.confname for x in self.head3list]
        for name in self.confnames:
            if len(name) != 14:
                print("      ERROR: %s is not a conformer name.")
                sys.exit()
            occurrence = self.confnames.count(name)
            if occurrence > 1:
                print("      ERROR: Conformer %s occurred %d times" % occurrence)
                sys.exit()
        return

env = Env()

if __name__ == "__main__":
    print("This is pymcce module.")
    print("Use pymonte.py to run mcce step 4.")