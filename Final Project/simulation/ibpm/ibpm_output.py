import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import savemat
data_path = './von_karman/'

filenames = sorted([os.path.join(data_path,f) for f in os.listdir(data_path) if f[-3:] == 'plt'])
timesteps = len(filenames)

U = np.zeros((449,199,timesteps))
V = np.zeros((449,199,timesteps))
W = np.zeros((449,199,timesteps))

start = time.time()

for timestep in range(timesteps):
    
    timestep_data = np.genfromtxt(filenames[timestep], delimiter=' ',skip_header=6)
    
    for i in range(449):
        
        for j in range(199):
            
            U[i,j,timestep] = timestep_data[i+449*j, 2]
            V[i,j,timestep] = timestep_data[i+449*j, 3]
            W[i,j,timestep] = timestep_data[i+449*j, 4]
            
    print('\rTimestep', timestep+1, 'of', timesteps, 'eta:', \
          int((timesteps-timestep-1)*(timestep+1)/(time.time()-start)), end = 's')
    

savemat('output.mat', {'U': U, 'V' : V, 'W': W})

# plot the data
plt.figure(figsize=(15, 8))

xx, yy = np.meshgrid(np.arange(449), np.arange(199))

for j in range(4):
    plt.subplot(2, 2, j+1)
    plt.pcolor(xx, yy, W[:, :, 20*j].T, cmap='coolwarm', vmin=-4, vmax=4)
    plt.colorbar()  # Adds a colorbar to each subplot to indicate the scale.

plt.show()