#!/usr/bin/python3

import MDAnalysis
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

topfile = "5cg.prmtop"
trjfile = "5cg-1ns.nc"

reftop = "5cg-ff-100ns.pdb"

# Read in topology and trajectories
trj = MDAnalysis.Universe(topfile,trjfile)
print(trj,len(trj.trajectory))
nframes = len(trj.trajectory)
natoms = len(trj.atoms)
nres = len(trj.residues)
nseg = len(trj.segments)

print("Nframes",nframes,"natoms",natoms,"nresidue",nres,"nsegments",nseg)
# Read in PDB as reference
ref = MDAnalysis.Universe(reftop)

# Select heavy-atom backbone
selatm = 'name CA C O N'

bb = trj.select_atoms(selatm)
nresbb = len(bb.residues)
print(len(bb.atoms),len(bb.residues))

# Select heavy atoms in PDB reference
refbb = ref.select_atoms(selatm)
print(len(bb.atoms),len(bb.residues))

#get reference posistions (subtracting COM)
refpos = refbb.translate(-refbb.center_of_mass())

val = np.zeros(nframes)
heatmap = np.zeros((nframes,nresbb+2))
fr = 0

for frame in trj.trajectory:
#    val[fr] = rmsd(bb.positions, refbb.positions, center=True, superposition=True)
#    print(rmsd(bb.positions,refbb.positions)) # raw rmsd
    trj0 = bb.translate(-bb.center_of_mass()) # move atoms to origin
    # print(rmsd(bb.positions,refbb.positions)) # rmsd traslated to com
    rotmat, rd  = align.rotation_matrix(trj0.positions,refpos.positions) #align with moved atom
    trj0 = bb.rotate(rotmat)
    rdr = rmsd(trj0.positions,refbb.positions) # rmsd translated and rotated
    val[fr] = rd
    heatmap[fr][0] = fr
    for ires in range(nresbb):
        ressel = "resid " + str(ires+1)
        refbbres = refbb.select_atoms(ressel)
        trj0res = trj0.select_atoms(ressel)

        #print(fr,ires,len(refbbres.atoms),len(trj0res.atoms))
        #        exit(1)
        if(len(refbbres.atoms) > 0):
            heatmap[fr][ires+1] = rmsd(trj0res.positions,refbbres.positions)
        else :
            print("Zero atoms in res?",fr,ires)
    
            
    # print(rd,rdr)
    # print(fr, val[fr])
    heatmap[fr][-1] = np.average(heatmap[fr][1:-2])
    #print(fr,val[fr],heatmap[fr][-1])
    fr+=1 # next frame

# To do
#rmsd without center and position
#center_of_mass
#align

print(np.average(val), np.average(heatmap[:,-1]),np.std(val))
#np.savetxt("heatmap.txt",np.column_stack(heatmap),fmt="%g")
np.savetxt("heatmap.txt",heatmap,fmt="%g")
# RMS-fitting
#rmsd(bb.positions, refbb.positions, center=True, superposition=True)

# Select all PA residues
# tail_pa = universe.select_atoms("resname PA")

# Select all 15th carbon of PA residues
# c15 = tail_pa.select_atoms("name C15")


