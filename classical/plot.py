import os, sys
import argparse
import h5py as hf
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--out_dir", dest="out_dir", default="figures", type=str, help="Directory that saving all resutls figures")
parser.add_argument("--data_name", dest="data_name", default="results.h5", type=str, help="Name of input data")
FLAGS = parser.parse_args()

if not os.path.exists(FLAGS.out_dir):
    os.makedirs(FLAGS.out_dir)

data = hf.File(FLAGS.data_name, "r")

T = data["T"][:]
E = data["E"][:] 
E2 = data["E2"][:] 
E4 = data["E4"][:] 
M = data["M"][:]
#absM = data["absM"][:]
M2 = data["M2"][:]
M4 = data["M4"][:]

Cv = (E2 - E**2) / T ** 2
X  = (M2 - M**2) / T 
R2 = M4 / M2 ** 2
# Binder Cumulant
U2 = 1.5 - 0.5 * R2


def plot_and_save_figure(path, name, T, observable):
    plt.plot(T, observable)
    plt.xlabel("T")
    plt.ylabel(name)
    plt.savefig(os.path.join(path, name + ".png"))
    plt.clf()

# PLOTTING

plot_and_save_figure(FLAGS.out_dir, "Binder Cumulant", T, U2)
plot_and_save_figure(FLAGS.out_dir, "Heat Capacity", T, Cv)
plot_and_save_figure(FLAGS.out_dir, "Susceptibility", T, X)
plot_and_save_figure(FLAGS.out_dir, "Energy", T, E)
plot_and_save_figure(FLAGS.out_dir, "Magnetization", T, M)