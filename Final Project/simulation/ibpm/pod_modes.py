import numpy as np
from mpl_toolkits.mplot3d import axes3d
from PDE_FIND import *
import scipy.io as sio
import itertools
import matplotlib.pyplot as plt
import pysindy as ps

data = sio.loadmat('output.mat')
steps = 101
n = 449
m = 199
W = data['W'].reshape(n,m,steps)
U = data['U'].reshape(n,m,steps)
V = data['V'].reshape(n,m,steps)

dt = 0.2
dx = 0.02
dy = 0.02

xmin = 100
xmax = 425
ymin = 15
ymax = 185

W = W[xmin:xmax, ymin:ymax, :]
U = U[xmin:xmax, ymin:ymax, :]
V = V[xmin:xmax, ymin:ymax, :]

n,m,steps = W.shape

## POD analysis

W = W.reshape(n*m, steps)
U = U.reshape(n*m, steps)
V = V.reshape(n*m, steps)

#Subtracting the mean from each dataset
U_mean = U.mean(axis=1, keepdims=True)
V_mean = V.mean(axis=1, keepdims=True)
W_mean = W.mean(axis=1, keepdims=True)
U_demeaned = U - U_mean
V_demeaned = V - V_mean
W_demeaned = W - W_mean

[U_modes, U_singular_values, U_temp] = np.linalg.svd(U_demeaned, full_matrices = False)
[V_modes, V_singular_values, V_temp] = np.linalg.svd(V_demeaned, full_matrices= False)
[W_modes, W_singular_values, W_temp] = np.linalg.svd(W_demeaned, full_matrices= False)

num_modes = 5

fig, axes = plt.subplots(1,num_modes, figsize = (18,6), dpi = 100)

for i in range(num_modes):
    mode_shape = W_modes[:, i].reshape(n, m)
    contour = axes[i].contourf(mode_shape, levels=100, cmap='viridis')
    fig.colorbar(contour, ax=axes[i])
    axes[i].set_title(f'W Mode {i+1}')
    axes[i].set_xlabel('x')
    axes[i].set_ylabel('y')
    axes[i].axis('tight')  # Makes the plot more compact

plt.tight_layout()
plt.show()

#Ploting the Mode Shapes for u and v
"""
num_modes = 5

fig, axes = plt.subplots(2,num_modes, figsize = (18,6), dpi = 100)

for i in range(num_modes):
    mode_shape = U_modes[:, i].reshape(n, m)
    contour = axes[0, i].contourf(mode_shape, levels=100, cmap='viridis')
    fig.colorbar(contour, ax=axes[0, i])
    axes[0, i].set_title(f'U Mode {i+1}')
    axes[0, i].set_xlabel('x')
    axes[0, i].set_ylabel('y')
    axes[0, i].axis('tight')  # Makes the plot more compact

# Plot modes for V
for i in range(num_modes):
    mode_shape = V_modes[:, i].reshape(n, m)
    contour = axes[1, i].contourf(mode_shape, levels=100, cmap='viridis')
    fig.colorbar(contour, ax=axes[1, i])
    axes[1, i].set_title(f'V Mode {i+1}')
    axes[1, i].set_xlabel('x')
    axes[1, i].set_ylabel('y')
    axes[1, i].axis('tight')  # Makes the plot more compact

plt.tight_layout()
plt.show()
"""