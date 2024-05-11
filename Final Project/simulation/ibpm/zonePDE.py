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

#SVD
W = W.reshape(n*m, steps)
U = U.reshape(n*m, steps)
V = V.reshape(n*m, steps)

uw, sigmaw, vw = np.linalg.svd(W, full_matrices = False); vw = vw.T
uu, sigmau, vu = np.linalg.svd(U, full_matrices=False); vu = vu.T
uv, sigmav, vv = np.linalg.svd(V, full_matrices=False); vv = vv.T

#plt.semilogy(sigmaw)
#plt.semilogy(sigmau)
#plt.semilogy(sigmav)
#plt.show()

Wn = W.reshape(n,m,steps)
Un = U.reshape(n,m,steps)
Vn = V.reshape(n,m,steps)

#Defining zones
zones = {
    'Zone11': {'xmin': 100, 'xmax': 150, 'ymin': 15, 'ymax': 50},
    'Zone12': {'xmin': 100, 'xmax': 150, 'ymin': 50, 'ymax': 100},
    'Zone13': {'xmin': 100, 'xmax': 150, 'ymin': 100, 'ymax': 150},
    'Zone14': {'xmin': 100, 'xmax': 150, 'ymin': 150, 'ymax': 185},

    'Zone21': {'xmin': 150, 'xmax': 200, 'ymin': 15, 'ymax': 50},
    'Zone22': {'xmin': 150, 'xmax': 200, 'ymin': 50, 'ymax': 100},
    'Zone23': {'xmin': 150, 'xmax': 200, 'ymin': 100, 'ymax': 150},
    'Zone24': {'xmin': 150, 'xmax': 200, 'ymin': 150, 'ymax': 185},

    'Zone31': {'xmin': 200, 'xmax': 250, 'ymin': 15, 'ymax': 50},
    'Zone32': {'xmin': 200, 'xmax': 250, 'ymin': 50, 'ymax': 100},
    'Zone33': {'xmin': 200, 'xmax': 250, 'ymin': 100, 'ymax': 150},
    'Zone34': {'xmin': 200, 'xmax': 250, 'ymin': 150, 'ymax': 185},

    'Zone41': {'xmin': 250, 'xmax': 300, 'ymin': 15, 'ymax': 50},
    'Zone42': {'xmin': 250, 'xmax': 300, 'ymin': 50, 'ymax': 100},
    'Zone43': {'xmin': 250, 'xmax': 300, 'ymin': 100, 'ymax': 150},
    'Zone44': {'xmin': 250, 'xmax': 300, 'ymin': 150, 'ymax': 185},

    'Zone51': {'xmin': 300, 'xmax': 350, 'ymin': 15, 'ymax': 50},
    'Zone52': {'xmin': 300, 'xmax': 350, 'ymin': 50, 'ymax': 100},
    'Zone53': {'xmin': 300, 'xmax': 350, 'ymin': 100, 'ymax': 150},
    'Zone54': {'xmin': 300, 'xmax': 350, 'ymin': 150, 'ymax': 185},

    'Zone61': {'xmin': 350, 'xmax': 400, 'ymin': 15, 'ymax': 50},
    'Zone62': {'xmin': 350, 'xmax': 400, 'ymin': 50, 'ymax': 100},
    'Zone63': {'xmin': 350, 'xmax': 400, 'ymin': 100, 'ymax': 150},
    'Zone64': {'xmin': 350, 'xmax': 400, 'ymin': 150, 'ymax': 185},
}

#Sampling points in zones
np.random.seed(0)
num_xy_per_zone = 3000
num_t = 40
num_points = num_xy_per_zone * num_t
boundary = 5
boundary_x = 10
points = {}
count = 0

points_zone = {zone: [] for zone in zones}
for zone, bounds in zones.items():
    count = 0
    while count < num_xy_per_zone:
        x = np.random.choice(range(bounds['xmin'] + boundary_x, bounds['xmax'] - boundary_x))
        y = np.random.choice(range(bounds['ymin'] + boundary, bounds['ymax'] - boundary))
        for t in range(num_t):
            points_zone[zone].append([x, y, 2 * t + 12])
        count += 1

#Run PDEfind
results = {}

N = 2*boundary-1  # odd number of points to use in fitting
Nx = 2*boundary_x-1  # odd number of points to use in fitting
deg = 5 # degree of polynomial to use

for zone, pts in points_zone.items():
    # Create data structures for derivatives and function values
    num_points = len(pts)
    w = np.zeros((num_points,1))
    u = np.zeros((num_points,1))
    v = np.zeros((num_points,1))
    wt = np.zeros((num_points,1))
    wx = np.zeros((num_points,1))
    wy = np.zeros((num_points,1))
    wxx = np.zeros((num_points,1))
    wxy = np.zeros((num_points,1))
    wyy = np.zeros((num_points,1))
        # Populate these arrays as before
    for i, point in enumerate(pts):
        [x, y, t] = point

        # Adjust x and y to be indices within the cropped array
        x = x - xmin  # Adjust x to start from 0 of the cropped range
        y = y - ymin  # Adjust y to start from 0 of the cropped range

        w[i] = Wn[x, y, t]
        u[i] = Un[x, y, t]
        v[i] = Vn[x, y, t]

        wt[i] = PolyDiffPoint(Wn[x,y,t-(N-1)//2:t+(N+1)//2], np.arange(N)*dt, deg, 1)[0]
    
        x_diff = PolyDiffPoint(Wn[x-(Nx-1)//2:x+(Nx+1)//2,y,t], np.arange(Nx)*dx, deg, 2)
        y_diff = PolyDiffPoint(Wn[x,y-(N-1)//2:y+(N+1)//2,t], np.arange(N)*dy, deg, 2)
        wx[i] = x_diff[0]
        wy[i] = y_diff[0]
    
        x_diff_yp = PolyDiffPoint(Wn[x-(Nx-1)//2:x+(Nx+1)//2,y+1,t], np.arange(Nx)*dx, deg, 2)
        x_diff_ym = PolyDiffPoint(Wn[x-(Nx-1)//2:x+(Nx+1)//2,y-1,t], np.arange(Nx)*dx, deg, 2)
    
        wxx[i] = x_diff[1]
        wxy[i] = (x_diff_yp[0]-x_diff_ym[0])/(2*dy)
        wyy[i] = y_diff[1]
       

    # Prepare the dataset and fit the model as previously done
    X_data = np.hstack([w, u, v])
    X_ders = np.hstack([np.ones((num_points,1)), wx, wy, wxx, wxy, wyy])
    X_ders_descr = ['','w_{x}', 'w_{y}','w_{xx}','w_{xy}','w_{yy}']
    X, description = build_Theta(X_data, X_ders, X_ders_descr, 2, data_description = ['w','u','v'])

    lam = 10**-5 
    d_tol = 5
    c = TrainSTRidge(X, wt, lam, d_tol)
    pde_result = print_pde(c, description, ut='w_t')

