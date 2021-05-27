import numpy as np
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import os


# Extract relevant data

base = os.getcwd()
path = base + '/pickle_out.txt'


in_file = open(path, 'rb')
A = pickle.load(in_file)
in_file.close()

_, solutions = A[0], A[1]


# Processing data (symmetrization along symmetry plane)

T2_tot = []
for i in range(len(solutions)):
    t1, t2 = solutions[i]
    T2_tot.append(t2)

T2 = np.asarray(T2_tot)


T2_r = np.flip(T2, axis=0)[1:]
T2_c = np.concatenate((T2_r, T2))

sistems = len(solutions)

add_0 = np.asarray([np.ones(2*sistems-1)])
T2_c = np.concatenate((add_0, T2_c.T), axis=0).T


# Preparing data for 3d plot

x = np.linspace(-(sistems-1), sistems-1, 2*sistems-1)*2.17/(sistems-1)  # converting to physical units

steps = solutions[0][0].shape[0]

y = np.linspace(0, steps-1, steps+1)*0.0001  # converting to physical units
X, Y = np.meshgrid(x, y)
X = X.T
Y = Y.T
Z = T2_c.real


# Plot the surface.

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize axes.

ax.set_zlim(0, 1.0)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zlabel(r'$Tr[\rho^2] $')
ax.set_xlabel(r'$\log_{10}(s)$')
ax.set_ylabel("t [u.a.]")

# Add color bar.

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

