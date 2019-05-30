import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


def plot_traje(fname, winkel):
    """Plotte die Datei, die in fname gegeben ist in einer 3D-Darstellung mit
    dem Azimuth(?)-Winkel winkel."""
    print("Plot "+fname)
    Y = np.genfromtxt(fname)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.view_init(winkel, 30)
    p = ax.scatter3D(Y[1, :], Y[2, :], Y[3, :], c=Y[0, :], cmap='viridis')
    cbar = fig.colorbar(p)
    cbar.set_label(r'Zeitschritte $t_n$', rotation=270, labelpad=10)
    fig.savefig(fname.split('.')[0]+'.jpg', dpi=600, bbox_inches='tight')
    fig.clf()


def plot_energie(fname):
    print("Plot "+fname)
    df = pd.read_csv(fname, decimal='.',
                     delimiter=' ')
    df.energie = df.energie - df.energie[0]
    fig = plt.figure()
    plt.plot(df.zeit, -df.energie, 'k-')
    plt.xlabel('Zeit')
    plt.ylabel('Energieabweichung')
    plt.yscale('log')
    fig.savefig(fname.split('.')[0]+'.jpg', dpi=600, bbox_inches='tight')
    fig.clf()


def plot_drehimpuls(fname):
    print("Plot "+fname)
    df = pd.read_csv(fname, decimal='.',
                     delimiter=' ')
    df.Lx = (df.Lx - df.Lx[0]).abs()
    df.Ly = (df.Ly - df.Ly[0]).abs()
    df.Lz = (df.Lz - df.Lz[0]).abs()
    fig = plt.figure()
    plt.plot(df.zeit, df.Lx, '-', label=r'$L_\mathrm{x}$')
    plt.plot(df.zeit, df.Ly, '-', label=r'$L_\mathrm{y}$')
    plt.plot(df.zeit, df.Lz, '-', label=r'$L_\mathrm{z}$')
    plt.xlabel('Zeit')
    plt.ylabel('Drehimpulsabweichung')
    plt.yscale('log')
    plt.legend(loc='best')
    fig.savefig(fname.split('.')[0]+'.jpg', dpi=600, bbox_inches='tight')
    fig.clf()


#  plot_traje('build/aufg2_a_ellipse.txt', 50)
#  plot_traje('build/aufg2_a_schmetterling.txt', 50)
plot_energie('build/aufg2_b_energie.txt')
plot_drehimpuls('build/aufg2_b_drehimpuls.txt')
#  plot_traje('build/aufg2_c_11.txt', 50)
#  plot_traje('build/aufg2_c_09.txt', 50)
