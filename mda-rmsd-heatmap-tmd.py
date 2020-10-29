#!/usr/bin/python3

import MDAnalysis
import sys
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

topfile = "5cg-test.prmtop"
trjfile = "5cg-test.nc"
reftop = "5cg-test.prmtop"
reftrj = "5cg-test.nc"

fres = 163 # first residue
lres = 394 # last residue

perresfile = "perres.txt"
heatmapfile = "heatmap.txt"

# Read in topology and trajectories
trj = MDAnalysis.Universe(topfile,trjfile)  # trajectory
ref = MDAnalysis.Universe(reftop,reftrj)    # reference

#print(trj,len(trj.trajectory))
nframes = len(trj.trajectory)
natoms = len(trj.atoms)
nres = len(trj.residues)
nseg = len(trj.segments)

print("Nframes:",nframes," natoms:",natoms,", nresidue:",nres,", nsegments:",nseg)

# Select heavy-atom backbone in trajectory
satom = 'name CA C O N and resid ' + str(fres) + ':' + str(lres)
bb = trj.select_atoms(satom) ###### get heavy-atom coordinates at first frame

print("Atoms Selected",len(bb.atoms)," Residue Selected",len(bb.residues))

# Select atoms for reference
refbb = ref.select_atoms(satom) #######
#lf = ref.trajectory[-1] # use last frame instead of first frame

# error check
if(len(bb.atoms) != len(refbb.atoms)): # do we have the same number of atoms?
   print("Error: # of atoms selected from trj do not match reference",len(bb.atoms),len(refbb.atoms))
   exit(1)
if(len(bb.residues) != len(refbb.residues)): # do we have the same number of residues
   print("Error: # of residues selected from trj do not match reference",len(bb.residues) != len(refbb.residues))
   exit(1)

#get reference posistions (subtracting COM)
ref0 = refbb.translate(-refbb.center_of_mass()) ####### move atoms to origin, renamed to ref0
ref0pos = ref0.positions ####### get ref0 coordinates
#print(refbb.center_of_mass(),ref0.center_of_mass())

# Allocate memory
val = np.zeros(nframes+1) # global rmsd
heatmap = np.zeros((nframes+1,nres+2)) # first column is frame number (+1 to offset python), last is average
fr = 1

for ires in range(nres+1):
    heatmap[0][ires]=ires # [0] is the first column, [ires] is the first row

for frame in trj.trajectory: # loop over configurations in trj
    trj0 = bb.translate(-bb.center_of_mass()) # move atoms to origin
    trj0pos = trj0.positions ####### assign trj0 positions to variable to use during frame iteration

    rotmat, rd  = align.rotation_matrix(trj0pos,ref0pos) # calculate rotation matrix for alignment
    trj0 = bb.rotate(rotmat) # do alignment
    rdr = rmsd(trj0pos,ref0pos) # global rmsd translated and rotated

    val[fr] = rd # global rmsd
    heatmap[fr][0] = fr # [fr is frame number] [0th row is residue number]
    for ires in range(fres, lres+1): # loop over residues
        ressel = "resid " + str(ires) # which atoms in residue ires
        refbbres = refbb.select_atoms(ressel) # select reference atoms in residue
        trj0res = trj0.select_atoms(ressel) # select trajectory atoms in trajectory

        if(len(refbbres.atoms) > 0): # error check that we have atoms!
            heatmap[fr][ires] = rmsd(trj0res.positions,refbbres.positions) # calculate the rmsd for this residue.
        else :
            print("Zero atoms in frame",fr," residue ",ires)
            
    heatmap[fr][-1] = np.average(heatmap[fr][1:-2]) # get average per res, should be close to global
    print(fr,val[fr],heatmap[fr][-1])
    fr+=1 # next frame

# create output

print(np.average(val), np.average(heatmap[:,-1]),np.std(val))

# header string
shead = sys.argv[0] + " " + topfile + " " + trjfile + " vs reference " + reftop + " " + reftrj 

print("Creating average per residue file: ",perresfile)
outfile = open(perresfile,'w')
outfile.write("# " + shead +"\n")
for i in range(nres):
    #print(i+1, np.average(heatmap[1:-1,i+1]))
    line = str(i+1) + ' ' + str(np.average(heatmap[1:-1,i+1])) + '\n'
    outfile.write(line)
outfile.close()

print("Creating heatmap file: ",heatmapfile)
shead = sys.argv[0] + " " + topfile + " " + trjfile + " vs reference " + reftop + " " + reftrj 
np.savetxt(heatmapfile,np.column_stack(heatmap),header=shead,fmt="%g")
#np.savetxt("heatmap-firstframe.txt",heatmap,fmt="%g")




