import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

datos = np.loadtxt("trayectoria_E_alpha.dat")

fig = figure()
ax = Axes3D(fig)
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.plot(datos[:,0], datos[:,1], datos[:, 2])
plt.savefig('trayectoria.pdf')
show()
