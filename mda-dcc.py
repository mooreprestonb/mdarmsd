#!/usr/bin/python3
# calculate dynamics cross correlation (DCC)
# deltaRi(t) = ri(t) - <ri(t)>
# DCC(ij) = <deltaRi(t) . deltaRj(t)>/[sqrt(<ri(t)^2>) sqrt(rj(t)^2)]

import MDAnalysis
import sys
import numpy as np

topfile = "5cg-test.prmtop"
trjfile = "5cg-test.nc"
report = 10

# 1-163 cytoplasmic domain
# 163-394 transmembrane domain

fres = 1 # first residue 
lres = 394 # last residue

# avgfile = "avg.txt"
dccmapfile = "dccmap4.txt"

# Read in topology and trajectories
trj = MDAnalysis.Universe(topfile,trjfile)  # trajectory

#print(trj,len(trj.trajectory))
nframes = len(trj.trajectory)
natoms = len(trj.atoms)
nres = len(trj.residues)
nseg = len(trj.segments)

print("Nframes:",nframes," natoms:",natoms,", nresidue:",nres,", nsegments:",nseg)

# Select CA backbone in trajectory
satom = 'name CA and resid ' + str(fres) + ':' + str(lres)
bb = trj.select_atoms(satom) ###### get CA coordinates

ncatms = len(bb.atoms)
print("Atoms Selected",len(bb.atoms)," Residue Selected",len(bb.residues))

#print(bb)
# Allocate memory
avgpos = np.zeros((ncatms,3)) # avgpos
dccmap = np.zeros((ncatms,ncatms)) # dcc map
ri2 = np.zeros(ncatms) # avgpos

fr = 1
#loop over frames to get avg pos
print("Getting average positions, report progress ever",report,"configs")
for frame in trj.trajectory: # loop over configurations in trj
   com = bb.center_of_mass()
   if(fr%report == 0):
      print("Frame of avg calc:",fr,com)
   ival = 0
   for atms in bb:
      avgpos[ival] += atms.position-com
      ival += 1
   fr+=1 # next frame
   
avgpos /= fr
#print(avgpos) # removed COM

print("Found average posistions from COM, looping over traj again to find dcc")
fr = 1  #loop over frames again to get dcc
print("Reporting progress every",report,"configs")
for frame in trj.trajectory: # loop over configurations in trj
   if(fr%report == 0):
      print("Frame of dcc calc:",fr,com)
   com = bb.center_of_mass()
   for ival in range(len(bb)):
      ri = bb[ival].position-com -avgpos[ival]
      ri2[ival] += np.dot(ri,ri)
      for jval in range(ival+1):
         rj = bb[jval].position-com - avgpos[jval]
         rd = np.dot(ri,rj)
         dccmap[ival][jval] += rd
#         dccmap[jval][ival] += rd
#         if(ival!=jval):
#            dccmap[jval][ival] += rd
   fr+=1 # next frame

# create output
print("Done:dcc... Normalizing and writing")
dccmap /= fr
r2 = np.sqrt(ri2/fr)

for ival in range(len(bb)):
   for jval in range(ival+1):
      rn = r2[ival]*r2[jval]
      dccmap[ival][jval] /= rn
      dccmap[jval][ival] = dccmap[ival][jval]  # make symmetric

print("Creating dccmap file: ",dccmapfile)
shead = sys.argv[0] + " " + topfile + " " + trjfile 
np.savetxt(dccmapfile,dccmap,header=shead,fmt="%g")
