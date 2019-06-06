import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


def plot(fname):
    Y = np.genfromtxt(fname)
    plt.plot(Y[:, 0], Y[:, 1])
    plt.xlabel(r'$v_0 / \mathrm{m s^{-1}}$')
    plt.ylabel(r'Betrag der vertikalen Geschwindigkeit / $\mathrm{m s^{-1}}$')
    plt.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    plt.clf()


plot("build/aufg2_a.txt")
