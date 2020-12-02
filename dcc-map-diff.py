#!/usr/bin/python3

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

mapfile1 = sys.argv[1]
mapfile2 = sys.argv[2]
#mapfile1 = sys.argv[1]
#mapfile2 = sys.argv[2]
hm1 = np.loadtxt(mapfile1) # load heatmap1 into matrix
hm2 = np.loadtxt(mapfile2) # load heatmap2 into matrix

hm = hm1-hm2
fig, ax = plt.subplots() # allocate plots
im = ax.imshow(hm, origin='lower') # 

# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax)

# Set title
ax.set_title("DCC")

#ax.set_yticks(np.arange(14))

plt.xlabel('Residue Number')
plt.ylabel('Residue Number')
plt.show()
