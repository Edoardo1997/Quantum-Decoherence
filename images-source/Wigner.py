import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
 

def gaussian(mu, sigma, x):
    return np.exp(-(x-mu)**2/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)

dx = 1
delta = 0.6
precision = 1000  # multiple of 50
limit = 3.2

delta_t = 0.7

x = np.linspace(-limit, limit, precision)
p = np.linspace(-limit, limit, precision)
X, P = np.meshgrid(x, p)


def Wigner(x, p, delta, x0, p0):
    return 1./np.pi*np.exp(-(x-x0)**2/delta**2-(p-p0)**2*delta**2)


disp = 2.

Wp = Wigner(X, P, delta, disp, 0)
Wm = Wigner(X, P, delta, -disp, 0)
I = Wigner(X, P, delta, 0, 0) * np.cos(P * 2 * disp)


Z = Wp + Wm + I
Z_f = Wp + Wm + I * np.exp(-3)  # 3 decoherence times


colortuple = ('black', 'white')
colors = np.empty(X.shape, dtype=str)


renorm = precision/50
for a in range(precision):
    for b in range(precision):
        colors[a, b] = colortuple[(int(a/renorm) + int(b/renorm)) % len(colortuple)]


fig = plt.figure()
ax = fig.gca(projection='3d')


Max = max(map(max, Z))

# Plot the surface.
surf = ax.plot_surface(X, P, Z, cmap=cm.seismic,
                       linewidth=0, antialiased=False, vmin=-Max, vmax=Max)

# Customize the z axis.
ax.set_zlim(-0.6, 0.6)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_zlabel(r'$W(x,p)$')
ax.set_xlabel(r'x')
ax.set_ylabel("p")


# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)


# Final
fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')


# Plot the surface.
surf = ax1.plot_surface(X, P, Z_f, cmap=cm.seismic,
                       linewidth=0, antialiased=False, vmin=-Max, vmax=Max)

# Customize the z axis.
ax1.set_zlim(-0.6, 0.6)
ax1.zaxis.set_major_locator(LinearLocator(10))
ax1.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax1.set_zlabel(r'$W(x,p)$')
ax1.set_xlabel(r'x')
ax1.set_ylabel("p")


# Add a color bar which maps values to colors.
fig1.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
