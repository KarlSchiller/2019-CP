import matplotlib.pyplot as plt
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df = pd.read_csv('build/aufg_a.txt', skiprows=1, decimal='.', delimiter=';')
# print(df.head())

plt.plot(df.x, df.pot, 'r.', label="Potential außerhalb des Würfels")
plt.legend(loc='best')
plt.xlabel(r'$x^*$')
plt.ylabel(r'$\Phi(x^*)$')
plt.savefig("build/aufg_a.pdf")
