import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(4.8, 4.8))  # default 6.4, 4.8


# Einlesen der Daten
#  df_10 = pd.read_csv('build/bild_10.txt', decimal='.', delimiter=';')
#  df_10 = df_10.loc[:, ~df_10.columns.str.contains('^Unnamed')]
# Plotten als Heatmaps
#  fig, ax = plt.subplots()
#  ax1 = sns.heatmap(df_10, linewidth=0.5, cmap='Greys')
#  plt.savefig('build/10.pdf')
#  plt.clf()

df_energy = pd.read_csv('build/aufg2-energy.txt', decimal='.', delimiter=';')
fig = plt.figure()
for i in range(0, 10):
    plt.plot(df_energy.n, df_energy['{}'.format(i)][0]-df_energy['{}'.format(i)], 'x', label='en {}'.format(i))
#  plt.plot(df_energy.n, df_energy['1'][0]-df_energy['1'], 'x', label='en 1')
plt.xlabel(r'$N$')
plt.yscale('log')
plt.legend(loc='best')
fig.savefig('build/aufg2-vergleich.pdf', bbox_inches='tight')
fig.clf()


def plot_eigenvalues():
    # plotte das Eigenwertspektrum
    df_eigen = pd.read_csv('build/aufg2-eigenvalues.txt', decimal='.', delimiter=';')

    indices = [i for i in range(360)]
    fig = plt.figure(figsize=(6.4, 4.8))
    plt.bar(indices, df_eigen.evalues, color='k')
    plt.xlabel('Eigenwert Nummer')
    plt.ylabel('Eigenwert')
    plt.yscale('log')
    fig.savefig('build/aufg2-eigenvalues.pdf', bbox_inches='tight')
    fig.clf()


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
