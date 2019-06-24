import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


def plot(fname, num_werte):
    """plotte die Magnetisierung mit Vergleichsplot."""
    print("Plot "+fname)
    Y = np.genfromtxt(fname)
    print(Y.shape)
    Y = Y.T
    print(Y.shape)
    magnet = np.array([])
    for i in range(num_werte):
        lol = Y[:, i]
        spinup = len(lol[lol > 0])
        # print(spinup)
        spindown = len(lol[lol < 0])
        # print(spinup - spindown)
        magnet = np.append(magnet, (spinup - spindown))

    maximum = max(magnet)
    magnet /= maximum
    h = np.linspace(-5, 5, num_werte)
    plt.plot(h, magnet, label='Numerisch')
    plt.plot(h, np.tanh(h), label='Analytisch')
    plt.xlabel(r'Äußeres Magnetfeld $H$')
    plt.ylabel(r'Magnetisierung $m$')
    plt.legend(loc='best')
    plt.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    plt.clf()


plot('build/aufg2.txt', 100)
