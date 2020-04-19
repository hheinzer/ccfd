#!/usr/bin/env python

import sys
import os
from glob import glob
import numpy as np

# check if there are two command line arguments
if len(sys.argv) < 2:
    print("ERROR: Not enough arguments!")
    sys.exit()

# get executable name and equation system
exe = sys.argv[1]
eqn = sys.argv[2]
print("Checking '%s', compiled with %s equation system"%(exe, eqn))

# get all ini files in the directory, sorted alphabetically
inis = sorted(glob("*.ini"))

# run the executable and compare the output with the target output
for ini in inis:
    print(" - %-18s: "%(ini), end="")
    os.system("%s %s > %s.log 2>/dev/null"%(exe, ini, ini.rpartition('.')[0]))

    # read target output
    targetName = sorted(glob("targetOutput/%s_*"%(ini.rpartition('.')[0])))[-1]
    targetDat = np.genfromtxt(targetName, skip_header=1, delimiter=',')[:,0:3]

    # read actual output
    outputName = sorted(glob("%s_*"%(ini.rpartition('.')[0])))[-1]
    outputDat = np.genfromtxt(outputName, skip_header=1, delimiter=',')[:,0:3]

    # calculate differences
    diff = np.sqrt(((targetDat - outputDat)**2).sum(axis=0))
    diff = np.sqrt(sum(diff**2))

    if diff == 0.0:
        print("good (%.4e)"%(diff))
    elif diff < 1e-5:
        print("ok   (%.4e)"%(diff))
    else:
        print("bad  (%.4e)"%(diff))
