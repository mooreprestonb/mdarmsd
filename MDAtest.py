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
    ALP = A.translate(-A.centor_of_mass())
    rmsdrot = align.rotation_matrix(ALP.posistion,AFPC)[1]
    rotn = A.translate(center_of_mass())
    rotn 
    print (ts.frame,rms.rmsd(ALP,AFP),rms.rmsd(ALPC,AFPC),rmsdrot,rotn)


