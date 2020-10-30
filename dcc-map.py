#!/usr/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

hm = np.loadtxt('dccmap.txt') # load heatmap into matrix

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



