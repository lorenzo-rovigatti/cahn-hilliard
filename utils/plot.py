from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
import numpy as np
import sys

if len(sys.argv) < 2:
    print("Usage is file1 [file2] [file3] [...]", file=sys.stderr)
    exit(1)

with open(sys.argv[1]) as input:
    N = len(input.readline().split())
    
for i, file in enumerate(sys.argv[1:]):
    plt.figure(i)
    plt.title(file)
    data = np.loadtxt(file)
    im = plt.imshow(data)
    plt.colorbar(im)
plt.show()
