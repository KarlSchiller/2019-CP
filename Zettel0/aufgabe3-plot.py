import matplotlib.pyplot as plt
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df = pd.read_csv('euler_vergleich.txt', skiprows=1, decimal='.', delimiter=';')
print(df.head())

plt.plot(df.t, df.normal, label=r"Normaler Euler")
plt.plot(df.t, df.symm, label=r"Symmetrischer Euler")
plt.plot(df.t, df.ana, label=r"$\exp(-t)$")
plt.xlabel(r"$t$")
plt.ylabel(r"$y(t)$")
plt.legend(loc='best')
plt.show()
