#!/usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

residues = [] # y-axis

for resnum in range(394):
    resnum += 1
    residues.append(resnum)
#print(residues)

frames = [] # x-axis

for frame in range(10000):
    frame += 1
    frames.append(frame)
#print(frames)

hm = np.loadtxt('heatmap-5cg-tmd-100ns.txt')

#hm = open('heatmap-5cg-100ns.txt','r')
#print(hm.readlines())

#movement = np.array(hm)
#print(movement)

fig, ax = plt.subplots()
im = ax.imshow(hm[1:,2:-2:20], origin='lower') # [remove 1st column (has labels), remove last row (has average)]

# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax)

# Set title
ax.set_title("Per-residue Receptor RMSD (5cg vs 5vai, 100ns)")

#ax.set_yticks(np.arange(14))
#ax.set_yticklabels(residues[::30])

plt.xlabel('Frames (x20)')
plt.ylabel('Residue Number')
plt.show()



