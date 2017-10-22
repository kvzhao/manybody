import os, sys

from os.path import isfile, join
from os import listdir

import argparse
import h5py as hf
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--data_dir", dest="data_dir", default="from_boromir", type=str, help="Folder of saved resutls")
parser.add_argument("--out_name", dest="out_name", default="out", type=str, help="Save figure name")
parser.add_argument("--alpha", dest="alpha", default=0.0, type=float, help="Universal constant alpha")
parser.add_argument("--beta", dest="beta", default=0.125, type=float, help="Universal constant beta")
parser.add_argument("--gamma", dest="gamma", default=1.75, type=float, help="Universal constant gamma")
parser.add_argument("--nu", dest="nu", default=1, type=float, help="Universal constant gamma")
parser.add_argument("--Tc", dest="Tc", default=2.27, type=float, help="Critical temperature of Ising Model")
parser.add_argument("--a", dest="a", default=1, type=float, help="A constant a")
parser.add_argument("--omega", dest="omega", default=0.001, type=float, help="A subleading exponent")
FLAGS = parser.parse_args()

# Read all data list from some folder
datafiles = [f for f in listdir(FLAGS.data_dir) if isfile(join(FLAGS.data_dir, f))]
datanames = []

def scaling_X(datafile):
    data = hf.File(join(FLAGS.data_dir, datafile), "r")
    L = int(df.rstrip(".h5").strip("L"))
    N = L ** 2
    T = data["T"][:]
    M = data["absM"][:]
    M2 = data["M2"][:]
    t = (T-FLAGS.Tc)/ FLAGS.Tc
    X = (M2 - M**2) / T * N
    xL = t * L ** (1/FLAGS.nu)
    yL = X * L ** (-FLAGS.gamma/FLAGS.nu) / (1 + FLAGS.a * L ** FLAGS.omega)
    return xL, yL

def scaling_Cv(datafile):
    data = hf.File(join(FLAGS.data_dir, datafile), "r")
    L = int(df.rstrip(".h5").strip("L"))
    N = L ** 2
    T = data["T"][:]
    E = data["E"][:]
    E2 = data["E2"][:]
    t = (T-FLAGS.Tc)/ FLAGS.Tc
    Cv = (E2 - E**2) / T ** 2 * N
    xL = t * L ** (1/FLAGS.nu)
    yL = Cv * L ** (-FLAGS.alpha/FLAGS.nu) #/ (1 + FLAGS.a * L ** FLAGS.omega)
    return xL, yL

def scaling_M(datafile):
    data = hf.File(join(FLAGS.data_dir, datafile), "r")
    L = int(df.rstrip(".h5").strip("L"))
    N = L ** 2
    T = data["T"][:]
    M = data["absM"][:]
    t = (T-FLAGS.Tc)/ FLAGS.Tc
    xL = t * L ** (1/FLAGS.nu)
    yL = M * L ** (FLAGS.beta/FLAGS.nu) #/ (1 + FLAGS.a * L ** FLAGS.omega)
    return xL, yL

for df in datafiles:
    xL, yL = scaling_X(df)
    datanames.append(df.rstrip(".h5"))
    plt.plot(xL, yL, '-*', linewidth=1.0)
titlename = " Susceptibility Scaling function - gamma={}, nu={}, Tc={}".format(FLAGS.gamma, FLAGS.nu, FLAGS.Tc)
plt.xlabel("xL")
plt.ylabel("yL")
plt.legend(datanames)
plt.grid(True)
plt.title(FLAGS.out_name + titlename)
plt.savefig(FLAGS.out_name + "_Susceptibility.png")
plt.clf()

for df in datafiles:
    xL, yL = scaling_Cv(df)
    datanames.append(df.rstrip(".h5"))
    plt.plot(xL, yL, '-*', linewidth=1.0)
titlename = " Heat Capacity Scaling function - alpha={}, nu={}, Tc={}".format(FLAGS.alpha, FLAGS.nu, FLAGS.Tc)
plt.xlabel("xL")
plt.ylabel("yL")
plt.legend(datanames)
plt.grid(True)
plt.title(FLAGS.out_name + titlename)
plt.savefig(FLAGS.out_name + "_HeatCapacity.png")
plt.clf()

for df in datafiles:
    xL, yL = scaling_M(df)
    datanames.append(df.rstrip(".h5"))
    plt.plot(xL, yL, '-*', linewidth=1.0)
titlename = " Magnetization Scaling function - beta={}, nu={}, Tc={}".format(FLAGS.beta, FLAGS.nu, FLAGS.Tc)
plt.xlabel("xL")
plt.ylabel("yL")
plt.legend(datanames)
plt.grid(True)
plt.title(FLAGS.out_name + titlename)
plt.savefig(FLAGS.out_name + "_Magnetization.png")