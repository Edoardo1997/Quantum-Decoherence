import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator


def gaussian(mu, sigma, x):
    return np.exp(-(x-mu)**2/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)


dx = 1
delta = 0.2
precision = 1000  # multiple of 50
limit = 2
alpha = 1.
delta_t = 0.7

x = np.linspace(-limit, limit, precision)
y = np.linspace(-limit, limit, precision)

gp = gaussian(dx, delta, x)
gm = gaussian(-dx, delta, x)
cat = (gp+gm)/np.sqrt(2)

X, Y = np.meshgrid(x, y)
Z = np.outer(cat, cat)
Z_f = Z * np.exp(-alpha * delta_t * (X-Y)**2)


# Create an empty array of strings with the same shape as the meshgrid
# and populate it with two colors in a checkerboard pattern.
colortuple = ('black', 'white')
colors = np.empty(X.shape, dtype=str)

renorm = precision/50
for a in range(precision):
    for b in range(precision):
        colors[a, b] = colortuple[(int(a/renorm) + int(b/renorm)) % len(colortuple)]

# initial plot
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, facecolors=colors, linewidth=0)


# final plot
fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
surf1 = ax1.plot_surface(X, Y, Z_f, facecolors=colors, linewidth=0)

# Customize the z axis.
ax.set_zlim(-2, 2)
ax.w_zaxis.set_major_locator(LinearLocator(6))

ax.set_title('Initial cat-state')
ax.set_zlabel("\u03C1"+"(x,y)")
ax.set_xlabel("x")
ax.set_ylabel("y")

ax1.set_zlim(-2, 2)
ax1.w_zaxis.set_major_locator(LinearLocator(6))
ax1.set_title('Final cat-state')
ax1.set_zlabel("\u03C1"+"(x,y)")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

plt.show()