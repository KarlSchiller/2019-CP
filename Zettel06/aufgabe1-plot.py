import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


def plot_newton():
    # plotte das Eigenwertspektrum
    df_newton = pd.read_csv('build/aufg1-newton.txt', decimal='.', delimiter=';')

    fig = plt.figure()
    #  plt.bar(indices, df_eigen.evalues, color='k')
    plt.plot(df_newton.i, df_newton.deriv, 'kx')
    plt.xlabel('Iteration')
    plt.ylabel(r"$f'\!\left(x\right)$")
    plt.yscale('log')
    fig.savefig('build/aufg1-newton.pdf', bbox_inches='tight')
    fig.clf()


def plot_bisection():
    # plotte das Eigenwertspektrum
    df_bisection = pd.read_csv('build/aufg1-bisection.txt', decimal='.', delimiter=';')

    fig = plt.figure()
    #  plt.bar(indices, df_eigen.evalues, color='k')
    plt.plot(df_bisection.i, df_bisection.upper - df_bisection.lower, 'kx')
    plt.xlabel('Iteration')
    plt.ylabel(r'$c - a$')
    plt.yscale('log')
    fig.savefig('build/aufg1-bisection.pdf', bbox_inches='tight')
    fig.clf()


plot_newton()
plot_bisection()
