import numpy as np
from PIL import Image
import sys

if len(sys.argv) != 4:
    print(f"Usage is {sys.argv[0]} conf_0 conf_1 conf_2", file=sys.stderr)
    exit(1)

with open(sys.argv[1]) as inp:
    N = len(inp.readline().split())

rho_max = 0
data = []
for i, f in enumerate(sys.argv[1:]):
    data.append(np.loadtxt(f))
    data_max = np.max(data[-1])
    if data_max > rho_max:
        rho_max = data_max

data = (np.dstack(data / rho_max) * 255.999).astype(np.uint8)
new_image = Image.fromarray(data)
new_image.save('saleh.png')
print("saleh.png printed", file=sys.stderr)
