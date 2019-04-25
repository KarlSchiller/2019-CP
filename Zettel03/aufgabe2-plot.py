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

df_energy = pd.read_csv('build/aufg2-energy.txt', decimal='.', delimiter=';')
fig = plt.figure()
#  n = np.linspace(10, 190, 19)
for i in range(0, 10):
    plt.plot(df_energy.n, df_energy['{}'.format(i)]-df_energy['{}'.format(i)][19], 'x', label='EW {}'.format(i))
    #  #  plt.plot(n, np.abs(np.diff(df_energy['{}'.format(i)])), 'x', label='en {}'.format(i))
#  plt.plot(df_energy.n, df_energy['1'][0]-df_energy['1'], 'x', label='en 1')
plt.xlabel(r'$N$')
plt.yscale('log')
plt.title('Energiedarstellung')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig('build/aufg2-energy.pdf', bbox_inches='tight')
fig.clf()

df_position = pd.read_csv('build/aufg2-position.txt', decimal='.', delimiter=';')
fig = plt.figure()
for i in range(0, 10):
    plt.plot(df_position.n, df_position['{}'.format(i)][19]-df_position['{}'.format(i)], 'x', label='EW {}'.format(i))
plt.xlabel(r'$N$')
plt.yscale('log')
plt.title('Ortsdarstellung')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
fig.savefig('build/aufg2-position.pdf', bbox_inches='tight')
fig.clf()
