#!/usr/bin/python3

import MDAnalysis
from MDAnalysis.analysis import align, rms
import numpy

# always start with a Universe, get topology and trajectory
u = MDAnalysis.Universe("gasby-interface.prmtop","gasby-5vai-100ns.nc")  

print(u)

A = u.select_atoms('backbone') # atoms of first frame
AFP = A.positions
AFPC = AFP - A.center_of_mass()

u.trajectory[-1] # select last frame

ALP = A.positions
ALPC = ALP - A.center_of_mass()

#print(AFP,ALP)
#raw rmsd
print(rms.rmsd(AFP,ALP))

#subtract center of mass
print(rms.rmsd(ALPC,AFPC))

print(align.rotation_matrix(ALPC,AFPC)) #[rotation matrix (ALPC to AFPC), rmsd]

print(A[0])
print(A[0].index,A[0].id,A[0].name,A[0].type)
print(A[0].resindex,A[0].resname,A[0].resid,A[0].resnum,A[0].segid)
print(A[0].position,A[0].mass,A[0].charge)
print(A[0].bonds[0])

#all frames
for ts in u.trajectory:
    AT = A.translate(-A.center_of_mass())
    ATrotmat, rmsdrot = align.rotation_matrix(AT.positions,AFPC)
    ATrot = AT.rotate(ATrotmat)
    print (ts.frame,rms.rmsd(AT.positions,AFPC),rms.rmsd(ATrot.positions,AFPC),rmsdrot)


