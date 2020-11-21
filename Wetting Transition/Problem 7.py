"""
Numerical solutions to d Omega/ d rho = 0 using Brent's method
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize

def func(rho, mu, epsilon, beta):
    return rho - (1 - rho) * np.exp(beta * (mu + 5 * epsilon * rho))

epsilon = 1
N = 1000
x = np.linspace(-5, 0, N) # mu / epsilon
y = np.linspace(0, 2, N) # 1 / (beta * epsilon)
roots = np.empty((len(x), len(y)))

for i in range(len(x)):
    for j in range(len(y)):
        roots[i, j] = scipy.optimize.brentq(func, 0, 1, args = (epsilon * x[i], epsilon, 1 / (y[j] * epsilon)))

fig, ax = plt.subplots()
plt.style.use("dark_background")
im = ax.imshow(roots, cmap = 'hot')

# Colorbar
clb = fig.colorbar(im)
clb.set_label('ρ', size = 20)

#ax.set_title('Heatmap of density')
ax.set_xticks(np.linspace(0, N, 6))
ax.set_yticks(np.linspace(0, N, 6))
ax.set_xticklabels(np.round_(np.linspace(0, 2, 6), 2), size = 15)
ax.set_yticklabels(np.round_(np.linspace(-5, 0, 6), 1), size = 15)
ax.set_xlabel('k T / ε', size = 20)
ax.set_ylabel('μ / ε', size = 20)

# Spinodal
xSpinodal = np.linspace(0, 0.6 * N, N)
ySpinodal = (N / 2) * np.ones(N)
plt.plot(xSpinodal, ySpinodal, 'c', linewidth = 3, label = 'BINODAL')

# Binodal
densities1 = 0.5 * (1 + np.sqrt(1 - ((4 * y) / 5)))
densities2 = 0.5 * (1 - np.sqrt(1 - ((4 * y) / 5)))
xSpinodal = np.linspace(0, N, N)
ySpinodal1 = (0.2 * N) * (y * np.log(densities1 / (1 - densities1)) - 5 * densities1) + N
ySpinodal2 = (0.2 * N) * (y * np.log(densities2 / (1 - densities2)) - 5 * densities2) + N
plt.plot(xSpinodal, ySpinodal1, 'b', linewidth = 3, label = 'SPINODAL')
plt.plot(xSpinodal, ySpinodal2, 'b', linewidth = 3)

plt.text(0.35 * N, 0.7 * N, 'LIQUID', color = 'k', size = 20)
plt.text(0.38 * N, 0.3 * N, 'GAS', color = 'w', size = 20)
plt.text(0 * N, 0.6 * N, 'UNSTABLE', color = 'w', size = 18)

plt.legend(prop={'size': 15})
fig.tight_layout()
plt.savefig("Density Heatmap.pdf")
plt.show()