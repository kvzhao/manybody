import os, sys

from os.path import isfile, join
from os import listdir

import argparse
import h5py as hf
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--observ", dest="observ", default="M2", type=str, help="Physical Observable")
parser.add_argument("--data_dir", dest="data_dir", default="from_boromir", type=str, help="Folder of saved resutls")
parser.add_argument("--out_name", dest="out_name", default="out", type=str, help="Save figure name")
FLAGS = parser.parse_args()

# Read all data list from some folder
datafiles = [f for f in listdir(FLAGS.data_dir) if isfile(join(FLAGS.data_dir, f))]
datanames = []

for df in datafiles:
    data = hf.File(join(FLAGS.data_dir, df), "r")
    T = data["T"][:]
    if FLAGS.observ == 'Cv':
        E = data["E"][:] 
        E2 = data["E2"][:] 
        observ = (E2 - E**2) / T ** 2
    elif FLAGS.observ == "X":
        #M = data["M"][:]
        M = data["absM"][:]
        M2 = data["M2"][:]
        observ  = (M2 - M**2) / T 
    elif FLAGS.observ == "U2":
        M2 = data["M2"][:]
        M4 = data["M4"][:]
        R2 = M4 / M2 ** 2
        # Binder Cumulant
        observ = 1.5 - 0.5 * R2
    else:
        observ = data[FLAGS.observ][:]
    datanames.append(df.rstrip(".h5"))
    plt.plot(T, observ, '-*', linewidth=1.0)

plt.xlabel("T")
plt.ylabel(FLAGS.observ)
plt.legend(datanames)
plt.grid(True)
plt.title(FLAGS.out_name)
plt.savefig(FLAGS.out_name + ".png")