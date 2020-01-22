#!/usr/bin/env python
import glob, subprocess, math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.transform import Rotation as R

#write out new pdbs after alignment
def make_pdb(array,name):
	with open(f'../PDB/aligned_{name[:-4]}','w') as f:
		with open(f'../PDB/{name}','r') as infile:
			i=0
			for line in infile:
				x = f'{array[i][0]:.3f}'
				y = f'{array[i][1]:.3f}'
				z = f'{array[i][2]:.3f}'
				f.write(f'{line[:27]}{x:>11}{y:>8}{z:>8} {line[55:]}')
				i+=1

#get the residue pair corresponding to the longest axis through
#the pocket and return it
def initialize (array):#, name):
	D=pdist(array)
	D=squareform(D)
	temp = np.where(D == D.max())
	return temp[0]


#find center of mass and translate coordinate array to be centered
#about the origin
def center(array):
	centered_array=np.zeros((array.shape[0],array.shape[1]))
	com=array.mean(0)
	for i in range(array.shape[0]):
		centered_array[i]=array[i]-com
	return centered_array


#align the longest pairwise vector with the z axis by doing rotations
#about 2d projections
def align(array,resid):
	aligned=np.zeros((array.shape[0],array.shape[1]))
	res1,res2=resid[0],resid[1]
	vector=array[res1] - array[res2]
	mag=np.linalg.norm(vector)
	norm=vector/mag

	xyproj=np.array([norm[0],norm[1],0])
	xzproj=np.array([norm[0],0,norm[2]])
	yzproj=np.array([0,norm[1],norm[2]])
	theta=np.arccos(np.dot(xyproj,[0,1,0]))
	psi=np.arccos(np.dot(xzproj,[1,0,0]))
	phi=np.arccos(np.dot(yzproj,[0,0,1]))
	rot_theta=R.from_euler('z',theta,degrees=False)
	rot_psi=R.from_euler('y',psi,degrees=False)
	rot_phi=R.from_euler('x',phi,degrees=False)

	al1=rot_theta.apply(array)
	al2=rot_psi.apply(al1)
	aligned=rot_phi.apply(al2)

	return aligned


#acquire names of pockets and coords of pocket in separate numpy
#matrices, name_array and coords respectively
name_array=[]
for name in glob.glob('../PDB/*_pock*'):
	name_array.append(name.split('/')[2])

for entry in name_array:
	coords = np.genfromtxt(f'../PDB/{entry}',dtype=float,usecols=(6,7,8))

#run all the functions
	centered=center(coords)
	pair=initialize(coords)
	aligned=align(centered,pair)
	make_pdb(aligned,entry)
