import matplotlib.pyplot as plt
import numpy as np
import os

if not os.path.isdir('build'):
    os.path.mkdir('build')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


def linear(x, m, n):
    return m*x+n

x = [0, 2.5, -6.3, 4, -3.2, 5.3, 10.1, 9.5, -5.4, 12.7]
y = [4, 4.3, -3.9, 6.5, 0.7, 8.6, 13, 9.9, -3.6, 15.1]

# aus aufgabe2.cpp
a = [0.959951, 2.65694]

x_plot = np.linspace(min(x), max(x), 1000)
plt.plot(x, y, 'kx', label='Messwerte')
plt.plot(x_plot, linear(x_plot, *a), 'b-', label='least squares')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='best')
plt.savefig('build/leastsquares.pdf')
plt.clf()
