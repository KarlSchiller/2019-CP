import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df_eigen = pd.read_csv('build/aufg3-eigenvalues.txt', decimal='.', delimiter=';')

indices = [i for i in range(360)]
plt.bar(indices, df_eigen.evalues, color='k')
plt.xlabel('Eigenwert Nummer')
plt.yscale('log')
plt.savefig('build/aufg3-eigenvalues.pdf')
plt.clf()
