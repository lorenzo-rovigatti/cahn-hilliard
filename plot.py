from matplotlib import pyplot as plt
import numpy as np

z = np.loadtxt("last.dat")
shw = plt.imshow(z)
plt.colorbar(shw)
plt.show()
