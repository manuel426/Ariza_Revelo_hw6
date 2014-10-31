#from show import show
import numpy as np
import matplotlib.pyplot as plt
import pylab

datos = np.loadtxt("trayectoria_E_alpha.dat")
x = datos[:,0]
y_1 = datos[:,1]

pylab.plot(x,y_1)
plt.show()
