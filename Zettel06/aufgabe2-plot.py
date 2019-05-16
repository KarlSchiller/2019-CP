import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits . mplot3d import Axes3D
import pandas as pd


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)


def plot_conjugate(x):

    df = pd.read_csv('build/conjugate.txt', skiprows=1, decimal='.',
                     delimiter=';')

    x1 = df.iloc[0::x, 0]
    y1 = df.iloc[0::x, 2]
    X, Y = np.meshgrid(x1, y1)
    werte = rosen(X, Y)
    plt.contour(X, Y, werte, cmap="winter")
    plt.colorbar()
    plt.plot(df.x1, df.x2, 'r.', markersize=1,
             label="Schritte")
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(np.min(x1)-0.005, np.max(x1)+0.005)
    plt.ylim(np.min(y1)-0.005, np.max(y1)+0.005)
    plt.savefig("build/conjugate.pdf")
    plt.clf()

    x = np.linspace(0, df.shape[0]-1, df.shape[0])
    plt.plot(x, df.ek, 'r.', label=r"$\epsilon_k$")
    plt.xlabel(r"Iterationen $k$")
    plt.ylabel(r"$\epsilon_k$")
    plt.legend(loc='best')
    plt.savefig('build/e_conjugate.pdf')


def plot_gradient(x):
    df = pd.read_csv('build/gradient.txt', skiprows=1, decimal='.',
                     delimiter=';')

    x1 = df.iloc[0::x, 0]
    y1 = df.iloc[0::x, 2]
    X, Y = np.meshgrid(x1, y1)
    werte = rosen(X, Y)
    plt.contour(X, Y, werte, cmap="winter")
    plt.colorbar()
    plt.plot(df.x1, df.x2, 'r.', markersize=1,
             label="Schritte")
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(-1.02, 1.02)
    plt.ylim(-1.02, 1.02)
    plt.savefig("build/gradient.pdf")
    plt.clf()

    x = np.linspace(0, df.shape[0]-1, df.shape[0])
    plt.plot(x, df.ek, 'r.', label=r"$\epsilon_k$")
    plt.xlabel(r"Iterationen $k$")
    plt.ylabel(r"$\epsilon_k$")
    plt.legend(loc='best')
    plt.savefig('build/e_gradient.pdf')


def plot_b1(x):
    df = pd.read_csv('build/b1.txt', skiprows=1, decimal='.',
                     delimiter=';')

    x1 = df.iloc[0::x, 0]
    y1 = df.iloc[0::x, 2]
    X, Y = np.meshgrid(x1, y1)
    werte = f(X, Y)
    plt.contour(X, Y, werte, cmap="winter")
    plt.colorbar()
    # ax.clabel(test, inline=1, fontsize=10)
    # plt.plot(x1, y1, 'r.', markersize=1,
    #         label="Schritte")
    plt.plot(1.5, 2.3, 'rx', label="Startpunkt", markersize=10)
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(np.min(x1)-0.005, np.max(x1)+0.005)
    plt.ylim(np.min(y1)-0.005, np.max(y1)+0.005)
    plt.savefig("build/b1.pdf")
    plt.clf()


def plot_b2(x):
    df = pd.read_csv('build/b2.txt', skiprows=1, decimal='.',
                     delimiter=';')

    x1 = df.iloc[0::x, 0]
    y1 = df.iloc[0::x, 2]
    X, Y = np.meshgrid(x1, y1)
    werte = f(X, Y)
    plt.contour(X, Y, werte, cmap="winter")
    plt.colorbar()
    # plt.plot(x1, y1, 'r.', markersize=1,
    #         label="Schritte")
    plt.plot(-1.7, -1.9, 'rx', label="Startpunkt", markersize=10)
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(np.min(x1)-0.005, np.max(x1)+0.005)
    plt.ylim(np.min(y1)-0.005, np.max(y1)+0.005)
    plt.savefig("build/b2.pdf")
    plt.clf()


def f(x1, x2):
    return 1/(1+np.exp(-10*(x1*x2-3)**2)/(x1**2 + x2**2))


def rosen(x1, x2):
    return (1-x1**2) + 100*(x2-x1**2)**2


# plot_conjugate(2)
plot_gradient(10)
plot_b1(100)
plot_b2(500)
# print(rosen(0.991117, 0.982278))
# print(rosen(1.00234, 0.133408))


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
