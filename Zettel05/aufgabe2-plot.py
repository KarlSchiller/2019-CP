import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


# Einlesen der Daten
#  df_10 = pd.read_csv('build/bild_10.txt', decimal='.', delimiter=';')
#  df_10 = df_10.loc[:, ~df_10.columns.str.contains('^Unnamed')]
# Plotten als Heatmaps
#  fig, ax = plt.subplots()
#  ax1 = sns.heatmap(df_10, linewidth=0.5, cmap='Greys')
#  plt.savefig('build/10.pdf')
#  plt.clf()


def plot_fourier():
    # plotte das Eigenwertspektrum
    df = pd.read_csv('build/test.txt', decimal='.', delimiter=';')

    plt.plot(df.x, df.real, 'kx')
    plt.savefig('build/test.pdf', bbox_inches='tight')
    plt.clf()
    #  indices = [i for i in range(360)]
    #  fig = plt.figure(figsize=(6.4, 4.8))
    #  plt.bar(indices, df_eigen.evalues, color='k')
    #  plt.xlabel('Eigenwert Nummer')
    #  plt.ylabel('Eigenwert')
    #  plt.yscale('log')
    #  fig.savefig('build/aufg2-eigenvalues.pdf', bbox_inches='tight')
    #  fig.clf()


def plot_testpicture():
    # plotte das Testbild und die Rekonstruktionen
    origpic, k200 = np.loadtxt('build/aufg2-k200.txt', delimiter='; ', unpack=True)
    origpic, k300 = np.loadtxt('build/aufg2-k300.txt', delimiter='; ', unpack=True)
    origpic = np.resize(np.array([origpic]), new_shape=(112, 92))
    k200 = np.resize(np.array([k200]), new_shape=(112, 92))
    k300 = np.resize(np.array([k300]), new_shape=(112, 92))

    plt.imshow(origpic, cmap='binary', vmin=0, vmax=255)
    plt.colorbar()
    plt.savefig('build/aufg2-testpic.pdf', bbox_inches='tight')
    plt.clf()

    plt.imshow(k200, cmap='binary', vmin=0, vmax=255)
    plt.colorbar()
    plt.savefig('build/aufg2-k200.pdf', bbox_inches='tight')
    plt.clf()

    plt.imshow(k300, cmap='binary', vmin=0, vmax=255)
    plt.colorbar()
    plt.savefig('build/aufg2-k300.pdf', bbox_inches='tight')
    plt.clf()


plot_fourier()
#  plot_testpicture()
