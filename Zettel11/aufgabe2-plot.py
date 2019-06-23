import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


def plot(fname):
    """plotte die Magnetisierung mit Vergleichsplot."""
    print("Plot "+fname)
    Y = np.genfromtxt(fname)
    print(Y.shape)
    Y = Y.T
    print(Y.shape)
    magnet = []
    for i in range(100):
        lol = Y[:, i]
        spinup = len(lol[lol > 0])
        print(spinup)
        spindown = len(lol[lol < 0])
        print(spindown)
        magnet.append((spinup - spindown)/100000)

    h = np.linspace(-5, 5, 100)
    plt.plot(h, magnet)
    plt.plot(h, np.tanh(h))
    plt.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    plt.clf()


plot('build/aufg2.txt')
