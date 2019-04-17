import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def plot_eigenvalues():
    # plotte das Eigenwertspektrum
    df_eigen = pd.read_csv('build/aufg2-eigenvalues.txt', decimal='.', delimiter=';')

    indices = [i for i in range(360)]
    plt.bar(indices, df_eigen.evalues, color='k')
    plt.xlabel('Eigenwert Nummer')
    plt.ylabel('Eigenwert')
    plt.yscale('log')
    plt.savefig('build/aufg2-eigenvalues.pdf')
    plt.clf()


def plot_testpicture():
    # plotte das Testbild und die Rekonstruktionen
    origpic, k200 = np.loadtxt('build/aufg2-k200.txt', delimiter='; ', unpack=True)
    origpic, k300 = np.loadtxt('build/aufg2-k300.txt', delimiter='; ', unpack=True)
    origpic = np.resize(np.array([origpic]), new_shape=(112, 92))
    k200 = np.resize(np.array([k200]), new_shape=(112, 92))
    k300 = np.resize(np.array([k300]), new_shape=(112, 92))

    plt.imshow(origpic, cmap='binary', vmin=0, vmax=255)
    plt.colorbar()
    plt.savefig('build/aufg2-testpic.pdf')
    plt.clf()

    plt.imshow(k200, cmap='binary', vmin=0, vmax=255)
    plt.colorbar()
    plt.savefig('build/aufg2-k200.pdf')
    plt.clf()

    plt.imshow(k300, cmap='binary', vmin=0, vmax=255)
    plt.colorbar()
    plt.savefig('build/aufg2-k300.pdf')
    plt.clf()


#  plot_eigenvalues()
plot_testpicture()
