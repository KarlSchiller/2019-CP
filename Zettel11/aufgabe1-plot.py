#  import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

tugreen = '#80BA26'
tuorange = '#E36913'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
#  plt.rc('figure', figsize=(4.8, 9))  # default 6.4, 4.8


def plot_config(fname):
    '''2D Plot aller Ortsvektoren und der Permutation'''
    df = pd.read_csv(fname, decimal='.', delimiter=' ')
    # sortiere Ortsvektoren gemaess gegebener permutation
    df.sort_values(by='perm', inplace=True)
    df.reset_index(inplace=True, drop=True)

    #  fig = plt.figure(figsize=(6.4, 4.8))
    fig = plt.figure()
    plt.plot(df.x, df.y, color=tugreen, marker='x', linestyle='', label='Ortsvektoren')
    plt.plot(df.x[0], df.y[0], color='r', marker='x', linestyle='', label='Start')
    plt.plot(df.x[-1:], df.y[-1:], color='b', marker='x', linestyle='', label='Ziel')
    plt.plot(df.x, df.y, color='k', linestyle='-', linewidth=0.5, label='Weg')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.legend(loc='best')
    fig.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    fig.clf()


plot_config('build/init.txt')
plot_config('build/orte.txt')

plot_config('build/d9S1e2.txt')
plot_config('build/d9S1e3.txt')
plot_config('build/d9S1e4.txt')

plot_config('build/d99S1e1.txt')
plot_config('build/d99S1e2.txt')
plot_config('build/d99S1e3.txt')

plot_config('build/d999S1e1.txt')
plot_config('build/d999S1e2.txt')

#  files = os.listdir('build')
#  files.sort()
#  for index,item in enumerate(files):
    #  if item[0] != 'd':
        #  files.pop(index)
    #  if item[-1] != 't':
        #  files.pop(index)
#  for name in files:
    #  plot_config('build/'+name)
