
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Einlesen der Daten
df_10 = pd.read_csv('build/bild_10.txt', decimal='.', delimiter=';')
df_10 = df_10.loc[:, ~df_10.columns.str.contains('^Unnamed')]
df_20 = pd.read_csv('build/bild_20.txt', decimal='.', delimiter=';')
df_20 = df_20.loc[:, ~df_20.columns.str.contains('^Unnamed')]
df_50 = pd.read_csv('build/bild_50.txt', decimal='.', delimiter=';')
df_50 = df_50.loc[:, ~df_50.columns.str.contains('^Unnamed')]

# Plotten als Heatmaps
fig, ax = plt.subplots()
ax1 = sns.heatmap(df_10, linewidth=0.5, cmap='Greys')
plt.savefig('build/10.pdf')
plt.clf()
ax2 = sns.heatmap(df_20, linewidth=0.5, cmap='Greys')
plt.savefig('build/20.pdf')
plt.clf()
ax3 = sns.heatmap(df_50, linewidth=0.5, cmap='Greys')
plt.savefig('build/50.pdf')
plt.clf()
