#!/usr/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

hm = np.loadtxt('heatmap-5cg-tmd-100ns.txt') # load heatmap into matrix

fig, ax = plt.subplots() # allocate plots
im = ax.imshow(hm[1:,2:-2:20], origin='lower') # [remove 1st column (has labels), remove last row (has average), stride 20]

# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax)

# Set title
ax.set_title("Per-residue Receptor RMSD (5cg vs 5vai, 100ns)")

#ax.set_yticks(np.arange(14))

plt.xlabel('Frames (x20)')
plt.ylabel('Residue Number')
plt.show()



