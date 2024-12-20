import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
# from heatFlux_osmosisAnnual.totalTemp.py import *

### 1D

x = np.arange(0, 10)
y = [10,9,14,13,10,12,9,16,14,12]
f = interpolate.interp1d(x, y)
print(x)
xnew = np.arange(0, 9, 0.1)
print(xnew)
ynew = f(xnew)   # use interpolation function returned by `interp1d`
print(ynew)
plt.plot(x, y, 'o', xnew, ynew, '-')
plt.savefig('interpolating_exponential_function.png'); plt.show()

exit()
### 2D

x = np.arange(-5.01, 5.01, 0.25) # time
y = np.arange(-5.01, 5.01, 0.25) # flux heat budget
xx, yy = np.meshgrid(x, y)

z = np.sin(xx**2+yy**2) ## gotm heat budget
f = interpolate.interp2d(x, y, z, kind='cubic')
#Now use the obtained interpolation function and plot the result:

xnew = np.arange(-5.01, 5.01, 1e-2)
ynew = np.arange(-5.01, 5.01, 1e-2)
znew = f(xnew, ynew)
print(znew)
plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
plt.show()
