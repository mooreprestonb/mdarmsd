#!/usr/bin/python3

import MDAnalysis
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

#topfile = "5cg.prmtop"
#trjfile = "5cg-1ns.nc"
# reftop = "5cg-ff-100ns.pdb"

# Read in topology and trajectories
trj = MDAnalysis.Universe('5vai-394-c1-glp1-ref.prmtop','receptor-5vai-c1-glp1-100ns.nc')
print(trj,len(trj.trajectory))
print(trj,len(trj.trajectory))
nframes = len(trj.trajectory)
natoms = len(trj.atoms)
nres = len(trj.residues)
fres = 163
lres = 394
nseg = len(trj.segments)

print("Nframes",nframes,"natoms",natoms,"nresidue",nres,"nsegments",nseg)
# Read in PDB as reference
# ref = MDAnalysis.Universe(reftop)

# Select heavy-atom backbone in trajectory
satom = 'name CA C O N and resid ' + str(fres) + ':' + str(lres)
bb = trj.select_atoms(satom) ###### get heavy-atom coordinates at first frame
bbpos = bb.positions ####### assign bb coordinates to variable to use later
nresbb = len(bb.residues)
#print(nresbb)
print(len(bb.atoms),len(bb.residues))

ref = MDAnalysis.Universe('../../5vai/5vai-394-ref.prmtop', '../../5vai/receptor-5vai-100ns.nc')
#lf = ref.trajectory[-1]



# Select heavy atoms in first frame reference
refbb = ref.select_atoms(satom) #######

#refbb.write("refbb.pdb")
#refbb = "refbb.pdb"
#refbb = MDAnalysis.Universe('refbb.pdb')
#refbb = refbb.select_atoms('name CA C O N')

#refbbpos = refbb.positions #######
print(len(bb.atoms),len(bb.residues))


#get reference posistions (subtracting COM)
print(refbb.center_of_mass())
ref0 = refbb.translate(-refbb.center_of_mass()) ####### move atoms to origin, renamed to ref0
ref0pos = ref0.positions ####### get ref0 coordinates
print(refbb.center_of_mass())


val = np.zeros(nframes+1)
heatmap = np.zeros((nframes+1,nres+2)) # first column is frame number (+1 to offset python), last is average
fr = 1
#print(heatmap)


for ires in range(nres+1):
    heatmap[0][ires]=ires # [0] is the first column, [ires] is the first row


for frame in trj.trajectory:
#    val[fr] = rmsd(bb.positions, refbb.positions, center=True, superposition=True)
#    print(rmsd(bb.positions,refbb.positions)) # raw rmsd
    trj0 = bb.translate(-bb.center_of_mass()) # move atoms to origin
    trj0pos = trj0.positions ####### assign trj0 positions to variable to use during frame iteration
    #print(frame)
    #print(trj.trajectory[fr])

    # print(rmsd(bb.positions,refbb.positions)) # rmsd traslated to com
    rotmat, rd  = align.rotation_matrix(trj0pos,ref0pos) #align with moved atom, ????? Why double assignment
    #print(rd)
    trj0 = bb.rotate(rotmat)
    rdr = rmsd(trj0pos,ref0pos) # rmsd translated and rotated, ####### ref0pos variable instead of refbb.postions
    #print(rdr)
    val[fr] = rd
    heatmap[fr][0] = fr # [fr is frame number] [0th row is residue number]
    for ires in range(fres, lres+1):
        ressel = "resid " + str(ires)
        #print(ressel)
        refbbres = refbb.select_atoms(ressel)
        #print(trj.trajectory[0,fr])
        #refbbres2 = ref0pos.select_atoms(ressel)
        #refbbres2 = refbbres.select_atoms(ressel)
        #print(refbbres)
        trj0res = trj0.select_atoms(ressel)
        #print(fr,ires,len(refbbres.atoms),len(trj0res.atoms))
        #        exit(1)
        if(len(refbbres.atoms) > 0):
            heatmap[fr][ires] = rmsd(trj0res.positions,refbbres.positions)
        else :
            print("Zero atoms in res?",fr,ires)
    #print(fr, val[fr])
    heatmap[fr][-1] = np.average(heatmap[fr][1:-2]) #?????? Why 2 subs after heatmap
    print(fr,val[fr],heatmap[fr][-1])
    fr+=1 # next frame
   #if fr==100:
     #   break

# TO-DO ITEM (COMPLETED)
#rmsd without center and position
#center_of_mass
#align
print(np.average(val), np.average(heatmap[:,-1]),np.std(val))
np.savetxt("heatmap-5cg-tmd-100ns.txt",np.column_stack(heatmap),fmt="%g")

#print(heatmap[1:99,5])

#outfile = open('perres-5cg-all-100ns.txt','w')
#for i in range(nresbb):
#    print(i+1, np.average(heatmap[1:-1,i+1]))
#    line = str(i+1) + ' ' + str(np.average(heatmap[1:-1,i+1])) + '\n'
#    outfile.write(line)
#outfile.close()

#np.savetxt("heatmap-firstframe.txt",heatmap,fmt="%g")
# RMS-fitting
#rmsd(bb.positions, refbb.positions, center=True, superposition=True)

# RMS fitting of whole trajectory
#alignment = align.AlignTraj(trj, ref, filename='rmsfit.dcd')
#alignment.run()

# Select all PA residues
# tail_pa = universe.select_atoms("resname PA")

# Select all 15th carbon of PA residues
# c15 = tail_pa.select_atoms("name C15")

# TO-DO ITEM
# Check hm values with VMD
# Modify code to use trajectory frame as oppose to loading PDB (COMPLETED)
# Use last frame as reference instead of first (COMPLETED)


