import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits import mplot3d

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)

Y = np.genfromtxt("build/ergebnis.txt")

fig = plt.figure()
ax = plt.axes(projection='3d')
p = ax.scatter3D(Y[1, :], Y[2, :], Y[3, :], c=Y[0, :], cmap='viridis')
fig.colorbar(p)
plt.savefig("build/ergebnis.pdf")
