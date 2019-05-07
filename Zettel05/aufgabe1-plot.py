import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)

# Einlesen der Datem
df = pd.read_csv('build/aufg_a.txt', skiprows=1, decimal='.', delimiter=';')
df_b = pd.read_csv('build/aufg_b.txt', skiprows=1, decimal='.', delimiter=';')
df_c_1 = pd.read_csv('build/aufg_c_1.txt', skiprows=1,
                     decimal='.', delimiter=';')
df_c_2 = pd.read_csv('build/aufg_c_2.txt', skiprows=1,
                     decimal='.', delimiter=';')


def fun(x):
    return 1/x


# Plotten der Daten
x = np.linspace(1, 8, 1000)
plt.plot(df.x, df.pot, 'r.', label="Potential außerhalb des Würfels")
plt.plot(x, fun(x), label=r'$\frac{1}{x}$')
plt.legend(loc='best')
plt.yscale('log')
plt.xlabel(r'$x^*$')
plt.ylabel(r'$\Phi(x^*)$')
plt.savefig("build/aufg_a.pdf")
plt.clf()

plt.plot(df_b.x, df_b.pot, 'r.', label="Potential innerhalb des Würfels")
plt.legend(loc='best')
plt.xlabel(r'$x^*$')
plt.ylabel(r'$\Phi(x^*)$')
plt.savefig("build/aufg_b.pdf")
plt.clf()

plt.plot(df_c_1.x, df_c_1.pot, 'r.', label="Potential außerhalb des Würfels")
plt.plot(x, fun(x), label=r'$\frac{1}{x}$')
plt.legend(loc='best')
plt.yscale('log')
plt.xlabel(r'$x^*$')
plt.ylabel(r'$\Phi(x^*)$')
plt.savefig("build/aufg_c_1.pdf")
plt.clf()

plt.plot(df_c_2.x, df_c_2.pot, 'r.', label="Potential innerhalb des Würfels")
plt.legend(loc='best')
plt.xlabel(r'$x^*$')
plt.ylabel(r'$\Phi(x^*)$')
plt.savefig("build/aufg_c_2.pdf")
