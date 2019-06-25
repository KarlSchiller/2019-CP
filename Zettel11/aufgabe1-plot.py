import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

tugreen = '#80BA26'
tuorange = '#E36913'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
#  plt.rc('figure', figsize=(4.8, 9))  # default 6.4, 4.8


def plot_config(fname):
    '''2D Plot aller Ortsvektoren und der Permutation'''
    df = pd.read_csv(fname, decimal='.', delimiter=' ')
    #  fig = plt.figure(figsize=(6.4, 4.8))
    fig = plt.figure()
    plt.plot(df.x, df.y, color=tugreen, marker='x', linestyle='', label='Ortsvektoren')
    #  plt.plot(df.x, df.y, color='k', linestyle='-', linewidth=0.5, label='Weg')
    #  TODO: Plot des Weges einbauen mit gegebener Permutation
    #  plt.xlim((0, 8))
    plt.xlabel(r'$x$')
    #  plt.ylim((0, 8))
    plt.ylabel(r'$y$')
    plt.legend(loc='best')
    fig.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    fig.clf()


plot_config('build/init.txt')
