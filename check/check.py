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
print("Checking ")

# get all ini files in the directory, sorted alphabetically
inis = sorted(glob("*.ini"))

# run the executable and compare the output with the target output
for ini in inis:
    print(" - %-15s:"%(ini))
    os.system("../%s %s > %s.log"%(exe, ini, ini.rpartition('.')[0]))
