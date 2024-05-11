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

# Sample a collection of data points, stay away from edges so I can just use centered finite differences.
np.random.seed(0)

num_xy = 2000
num_t = 40
num_points = num_xy * num_t
boundary = 5
boundary_x = 10
points = {}
count = 0

for p in range(num_xy):
    x = np.random.choice(np.arange(boundary_x,n-boundary_x),1)[0]
    y = np.random.choice(np.arange(boundary,m-boundary),1)[0]
    for t in range(num_t):
        points[count] = [x,y,2*t+12]
        count = count + 1

# Take up to second order derivatives.
w = np.zeros((num_points,1))
u = np.zeros((num_points,1))
v = np.zeros((num_points,1))
wt = np.zeros((num_points,1))
wx = np.zeros((num_points,1))
wy = np.zeros((num_points,1))
wxx = np.zeros((num_points,1))
wxy = np.zeros((num_points,1))
wyy = np.zeros((num_points,1))

N = 2*boundary-1  # odd number of points to use in fitting
Nx = 2*boundary_x-1  # odd number of points to use in fitting
deg = 5 # degree of polynomial to use

'''
for p in points.keys():
    
    [x,y,t] = points[p]
    w[p] = Wn[x,y,t]
    u[p] = Un[x,y,t]
    v[p] = Vn[x,y,t]
    
    # Calculate the sizes of the inputs
    size_Wn = len(Wn[x,y,t-(N-1)//2:t+(N+1)//2])
    size_time = len(np.arange(N)*dt)
    
    # Print the sizes
    print("Size of Wn slice: ", size_Wn)
    print("Size of time array: ", size_time)
'''

for p in points.keys():
    
    [x,y,t] = points[p]
    w[p] = Wn[x,y,t]
    u[p] = Un[x,y,t]
    v[p] = Vn[x,y,t]
    
    wt[p] = PolyDiffPoint(Wn[x,y,t-(N-1)//2:t+(N+1)//2], np.arange(N)*dt, deg, 1)[0]
    
    x_diff = PolyDiffPoint(Wn[x-(Nx-1)//2:x+(Nx+1)//2,y,t], np.arange(Nx)*dx, deg, 2)
    y_diff = PolyDiffPoint(Wn[x,y-(N-1)//2:y+(N+1)//2,t], np.arange(N)*dy, deg, 2)
    wx[p] = x_diff[0]
    wy[p] = y_diff[0]
    
    x_diff_yp = PolyDiffPoint(Wn[x-(Nx-1)//2:x+(Nx+1)//2,y+1,t], np.arange(Nx)*dx, deg, 2)
    x_diff_ym = PolyDiffPoint(Wn[x-(Nx-1)//2:x+(Nx+1)//2,y-1,t], np.arange(Nx)*dx, deg, 2)
    
    wxx[p] = x_diff[1]
    wxy[p] = (x_diff_yp[0]-x_diff_ym[0])/(2*dy)
    wyy[p] = y_diff[1]

# Form a huge matrix using up to quadratic polynomials in all variables.
X_data = np.hstack([w,u,v])
X_ders = np.hstack([np.ones((num_points,1)), wx, wy, wxx, wxy, wyy])
X_ders_descr = ['','w_{x}', 'w_{y}','w_{xx}','w_{xy}','w_{yy}']
X, description = build_Theta(X_data, X_ders, X_ders_descr, 2, data_description = ['w','u','v'])
print('Candidate terms for PDE')
print(['1']+description[1:])


lam = 10**-5
d_tol = 5
c = TrainSTRidge(X,wt,lam,d_tol)
print_pde(c, description, ut = 'w_t')

