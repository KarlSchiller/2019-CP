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
    maximum = max(Y)
    Y /= maximum
    h = np.linspace(-5, 5, num_werte)
    plt.plot(h, Y, linewidth=5, label='Numerisch')
    plt.plot(h, np.tanh(h), label='Analytisch')
    plt.xlabel(r'Äußeres Magnetfeld $H$')
    plt.ylabel(r'Magnetisierung $m$')
    plt.legend(loc='best')
    plt.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    plt.clf()


plot('build/aufg2.txt', 10**4)
plot('build/aufg2_100.txt', 10**2)
