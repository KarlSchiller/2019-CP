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


def plot_particles(fname):
    '''2D Plot aller Teilchen'''
    df = pd.read_csv(fname, decimal='.', delimiter=' ')
    fig = plt.figure(figsize=(4.8, 4.8))
    plt.plot(df.x, df.y, 'kx', label='Teilchen')
    plt.xlim((0, 8))
    plt.xlabel(r'$x$')
    plt.ylim((0, 8))
    plt.ylabel(r'$y$')
    fig.savefig(fname.split('.')[0]+'.jpg', bbox_inches='tight')
    fig.clf()


def plot_equi(fname):
    '''Aequilibrierungsphase, plotte
    Schwerpunktsgeschwindigkeit
    kinetische Energie
    Temperatur'''
    df = pd.read_csv(fname, decimal='.', delimiter=' ')

    # Schwerpunktsgeschwindigkeit x-Komponente
    fig = plt.figure(figsize=(6.4, 4.8))
    plt.plot(df.t, df.vsx, color=tugreen, linestyle='-')
    plt.xlabel(r'$t$ ($\tau$)')
    plt.ylabel(r'$v_\mathrm{cms, x}$ ($\frac{\sigma}{\tau}$)')
    fig.savefig(fname.split('.')[0]+'_vSx.jpg', bbox_inches='tight')
    fig.clf()

    # Schwerpunktsgeschwindigkeit y-Komponente
    fig = plt.figure(figsize=(6.4, 4.8))
    plt.plot(df.t, df.vsy, color=tugreen, linestyle='-')
    plt.xlabel(r'$t$ ($\tau$)')
    plt.ylabel(r'$v_\mathrm{cms, y}$ ($\frac{\sigma}{\tau}$)')
    fig.savefig(fname.split('.')[0]+'_vSy.jpg', bbox_inches='tight')
    fig.clf()

    # Energie
    fig = plt.figure(figsize=(6.4, 4.8))
    plt.plot(df.t, df.Ekin, color=tugreen, linestyle='-', label=r'$E_\mathrm{kin}$')
    plt.plot(df.t, df.Epot, color=tuorange, linestyle='-', label=r'$E_\mathrm{pot}$')
    plt.plot(df.t, df.Ekin+df.Epot, color='r', linestyle='-', label=r'$E_\mathrm{ges}$')
    plt.xlabel(r'$t$ ($\tau$)')
    plt.ylabel(r'$E$ ($\varepsilon$)')
    plt.legend(loc='best')
    fig.savefig(fname.split('.')[0]+'_energy.jpg', bbox_inches='tight')
    fig.clf()

    # Temperatur
    fig = plt.figure(figsize=(6.4, 4.8))
    plt.plot(df.t, df.Temp, color=tugreen, linestyle='-')
    plt.xlabel(r'$t$ ($\tau$)')
    plt.ylabel(r'$T$ ($\frac{\varepsilon}{k_\mathrm{B}}$)')
    fig.savefig(fname.split('.')[0]+'_temp.jpg', bbox_inches='tight')
    fig.clf()


def plot_measurement(fname):
    '''Messungsphase, plotte
    Temperatur'''
    df = pd.read_csv(fname, decimal='.', delimiter=' ')
    print(fname)
    print('Temperaturmittel     ', df.Temp.mean(), ' +- ', df.Temp.std())
    #  Falsch verstanden: Hier wäre das Temperaturmittel über jeweils 1e4 Werte berechnet
    #  tmean = np.array([])
    #  Tmean = np.array([])
    #  for i in range(0, df.Temp.size, 10000):
        #  tmean = np.append(tmean, df.t[i:i+10000].mean())
        #  Tmean = np.append(Tmean, df.Temp[i:i+10000].mean())

    # Temperatur
    fig = plt.figure(figsize=(6.4, 4.8))
    #  plt.plot(tmean, Tmean, color=tugreen, linestyle='-')
    plt.plot(df.t, df.Temp, color=tugreen, linestyle='-')
    #  plt.xlabel(r'$t_\mathrm{mean}$ ($\tau$)')
    #  plt.ylabel(r'$T_\mathrm{mean}$ ($\frac{\varepsilon}{k_\mathrm{B}}$)')
    plt.xlabel(r'$t$ ($\tau$)')
    plt.ylabel(r'$T$ ($\frac{\varepsilon}{k_\mathrm{B}}$)')
    fig.savefig(fname.split('.')[0]+'_temp.jpg', bbox_inches='tight')
    fig.clf()


# Aufgabenteil a)
plot_particles('build/init.txt')  # Verifikation der Init-Methode
# Aufgabenteil b)
plot_equi('build/equilibration.txt')
plot_measurement('build/messung_T1.txt')
# Aufgabenteil c)
plot_measurement('build/messung_T1e2.txt')
plot_measurement('build/messung_T1e-2.txt')
# Aufgabenteil d)
plot_equi('build/thermostat_equi.txt')
plot_measurement('build/thermostat_mess.txt')
plot_particles('build/endkonfig.txt')  # Endkonfiguration
