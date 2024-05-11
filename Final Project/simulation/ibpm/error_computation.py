import numpy as np

err = abs(np.array([(0.009952-0.01)*100/0.01, (0.009870-0.01)*100/0.01, (-0.990560+1)*100, (-0.986379+1)*100]))

print("Error using PDE-FIND to identify Navier-Stokes:\n")
print("Mean parameter error:", np.mean(err), '%')
print("Standard deviation of parameter error:", np.std(err), '%')