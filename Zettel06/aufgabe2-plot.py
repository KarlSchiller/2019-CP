# import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits . mplot3d import Axes3D
import pandas as pd


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)


def plot_conjugate():

    df_conjugate = pd.read_csv('build/conjugate.txt', skiprows=1, decimal='.',
                               delimiter=';')

    plt.plot(df_conjugate.x1, df_conjugate.x2, 'r.', markersize=1,
             label="Schritte")
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.savefig("build/conjugate.pdf")
    plt.clf()


def plot_gradient():
    df_gradient = pd.read_csv('build/gradient.txt', skiprows=1, decimal='.',
                              delimiter=';')

    plt.plot(df_gradient.x1, df_gradient.x2, 'r.', markersize=1,
             label="Schritte")
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.savefig("build/gradient.pdf")
    plt.clf()


def plot_b1():
    df_gradient = pd.read_csv('build/b1.txt', skiprows=1, decimal='.',
                              delimiter=';')

    plt.plot(df_gradient.x1, df_gradient.x2, 'r.', markersize=1,
             label="Schritte")
    plt.plot(1.5, 2.3, 'bx', label="Startpunkt", markersize=10)
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.savefig("build/b1.pdf")
    plt.clf()


def plot_b2():
    df_gradient = pd.read_csv('build/b2.txt', skiprows=1, decimal='.',
                              delimiter=';')

    plt.plot(df_gradient.x1, df_gradient.x2, 'r.', markersize=1,
             label="Schritte")
    plt.plot(-1.7, -1.9, 'bx', label="Startpunkt", markersize=10)
    plt.legend(loc='best')
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.savefig("build/b2.pdf")
    plt.clf()


plot_conjugate()
plot_gradient()
plot_b1()
plot_b2()
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
