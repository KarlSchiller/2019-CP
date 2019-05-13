import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits . mplot3d import Axes3D


def f(x1, x2):
    return (1-x1**2) + 100*(x2-x1**2)**2


x1 = np.linspace(0, 2)
x2 = np.linspace(0, 2)
X1, X2 = np.meshgrid(x1, x2)
werte = f(X1, X2)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(
    X1, X2, werte,
    lw=0,
    s=5,
)
plt.savefig('build/test.pdf')
