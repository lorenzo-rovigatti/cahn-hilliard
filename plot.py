from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.animation as animation
import numpy as np
import sys

if len(sys.argv) < 3:
    print("Usage is N file", file=sys.stderr)
    exit(1)
    
N = int(sys.argv[1])
file = sys.argv[2]

fig, ax = plt.subplots()

frames = []
with open(file) as input_file:
    def list_chunks(lines, n):
        for i in range(0, len(lines), n):
            yield lines[i:i + n]
    
    for lines in list_chunks(input_file.readlines(), N):
        data = np.loadtxt(lines)
        frames.append(data)


# make sure that the color bar handles the values found at the end of the simulation
norm = colors.Normalize(vmin=frames[-1].min(), vmax=frames[-1].max())
image = plt.imshow(frames[0], norm=norm)
cbar = fig.colorbar(image, label="$\psi$")

def update_figure(j):
    image.set_data(frames[j])
    fig.canvas.flush_events()
    
    return [image]

anim = animation.FuncAnimation(fig, update_figure, frames=len(frames), interval=100, blit=True, repeat=False)
anim.save(file + ".mp4") # requires ffmpg
plt.show()
