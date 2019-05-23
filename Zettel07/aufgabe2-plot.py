import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.formatter', use_locale=True)
plt.rc('figure', figsize=(6.4, 4.8))  # default 6.4, 4.8


def plot_poly(fname):
    """Plotte Polygonzug, der in Datei @fname gespeichert ist"""
    df = pd.read_csv(fname, decimal='.', delimiter=' ', header=0)

    fig = plt.figure()
    plt.plot(df.x, df.y, 'kx')
    plt.xlabel('x')
    plt.ylabel('y')
    fig.savefig(fname.split('.')[0]+'.pdf', bbox_inches='tight')
    fig.clf()


plot_poly('build/aufg2-l1-start.txt')
plot_poly('build/aufg2-l2-start.txt')
