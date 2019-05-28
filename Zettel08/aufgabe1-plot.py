import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)


def plot_traje(fname, winkel):
    """Plotte die Datei, die in fname gegeben ist in einer 3D-Darstellung mit
    dem Azimuth(?)-Winkel winkel."""
    Y = np.genfromtxt(fname)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.view_init(winkel, 30)
    p = ax.scatter3D(Y[1, :], Y[2, :], Y[3, :], c=Y[0, :], cmap='viridis')
    cbar = fig.colorbar(p)
    cbar.set_label(r'Zeitschritte $t_n$', rotation=270, labelpad=10)
    fig.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    fig.clf()


plot_traje('build/a_harm.txt', 70)
plot_traje('build/a_unharm.txt', 70)
