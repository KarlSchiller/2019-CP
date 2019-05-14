# import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits . mplot3d import Axes3D
import pandas as pd


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)

df_conjugate = pd.read_csv('build/conjugate.txt', skiprows=1, decimal='.',
                           delimiter=';')

plt.plot(df_conjugate.x1, df_conjugate.x2, 'r.', markersize=1,
         label="Schritte")
plt.legend(loc='best')
plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")
plt.savefig("build/conjugate.pdf")
plt.clf()

df_gradient = pd.read_csv('build/gradient.txt', skiprows=1, decimal='.',
                          delimiter=';')

plt.plot(df_gradient.x1, df_gradient.x2, 'r.', markersize=1,
         label="Schritte")
plt.legend(loc='best')
plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")
plt.savefig("build/gradient.pdf")
# def f(x1, x2):
#     return (1-x1**2) + 100*(x2-x1**2)**2
#
#
# x1 = np.linspace(0, 2)
# x2 = np.linspace(0, 2)
# X1, X2 = np.meshgrid(x1, x2)
# werte = f(X1, X2)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(
#     X1, X2, werte,
#     lw=0,
#     s=5,
# )
# plt.savefig('build/test.pdf')
